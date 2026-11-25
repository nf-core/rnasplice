#!/usr/bin/env python

# Summarise the result tables of the splicing analyses into MultiQC custom-content files:
# General Statistics columns, and the sections of a Splicing parent section, whose text,
# headers and tool tables come from the section config.

import argparse
import html
import json
import os
import platform
import shlex
import sys
from collections import Counter, defaultdict, namedtuple
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import yaml

TAB = "\\t"
LABELS = json.loads(r"""${groovy.json.JsonOutput.toJson(labels)}""")
STAGED = json.loads(
    r"""${groovy.json.JsonOutput.toJson((results instanceof List ? results : [results]).collect { result -> result.toString() })}"""
)
EXT_ARGS = json.loads(r"""${groovy.json.JsonOutput.toJson([task.ext.args ?: ''])}""")[0]
PROCESS = "${task.process}"
with open("${section_config}") as handle:
    CONFIG = yaml.safe_load(handle)
SECTIONS, LIMITS, METHODS, UNITS = (
    CONFIG[key] for key in ("sections", "limits", "methods", "units")
)
ROUTES = ["STAR", "Salmon", ""]
DTU = {"DRIMSeq", "IsoformSwitchAnalyzeR"}
VARIANCE = {
    "DEXSeq DEU": "DEXSeq DEU",
    "edgeR DEU": "edgeR DEU",
    "DRIMSeq": "DEXSeq DTU",
    "SUPPA2": "SUPPA2",
}
GROUP = ["method", "route", "unit", "contrast"]
Contrast = namedtuple("Contrast", "treatment control")
Failure = namedtuple("Failure", "name filename error")


@dataclass
class Run:
    args: dict
    contrasts: dict = field(default_factory=dict)
    condition_of: dict = field(default_factory=dict)
    frames: list = field(default_factory=list)
    counts: list = field(default_factory=list)
    general_stats: dict = field(default_factory=lambda: defaultdict(dict))
    splices: dict = field(default_factory=dict)
    psi: dict = field(default_factory=dict)
    leafcutter: dict = field(default_factory=dict)
    failures: list = field(default_factory=list)
    ran: set = field(default_factory=set)

    def __post_init__(self):
        self.genes = Genes()

    def samples(self, condition):
        return {sample for sample, of in self.condition_of.items() if of == condition}

    def sizes(self):
        """Samples of the treatment and control of each contrast whose conditions are in the samplesheet."""
        present = set(self.condition_of.values())
        return {
            name: tuple(len(self.samples(c)) for c in pair)
            for name, pair in self.contrasts.items()
            if set(pair) <= present
        }


def parse_args(ext_args):
    parser = argparse.ArgumentParser()
    for name, default in CONFIG["args"].items():
        kind = (
            (lambda value: value.lower() == "true")
            if isinstance(default, bool)
            else type(default)
        )
        parser.add_argument(f"--{name}", type=kind, default=default)
    args = vars(parser.parse_args(shlex.split(ext_args)))
    return (
        args
        | LIMITS
        | {
            "suppa_correction": "gene-corrected"
            if args["suppa_gene_correction"]
            else "uncorrected"
        }
    )


def read_table(path, sep=TAB, **kwargs):
    """A table as strings, as written, so that numbers() treats every missing value alike."""
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, **kwargs)


def numbers(series):
    return pd.to_numeric(series, errors="coerce")


def features(method, unit, significant, counted=True, **columns):
    columns = {
        key: value
        if value is None or isinstance(value, str)
        else pd.Series(list(value), dtype=object).to_numpy()
        for key, value in columns.items()
    }
    return pd.DataFrame(
        {
            "method": method,
            "unit": unit,
            "significant": np.asarray(significant, dtype=bool),
            "counted": counted,
            **columns,
        }
    )


def with_route(name, route, template="{name} ({route})"):
    return template.format(name=name, route=CONFIG["routes"][route]) if route else name


def first_route(routes):
    """The STAR and Salmon routes quantify the same samples, so per-sample summaries show one."""
    return next((route for route in ROUTES if route in routes), None)


def message(section, key, **values):
    return SECTIONS[section]["messages"][key].format(**values)


def emit(section_id, plot_type, data, ctx=None, description=None, **extra):
    section, ctx = SECTIONS.get(section_id, {}), ctx or {}
    content = {"id": section_id}
    if section:
        references = (
            "\\n\\n" + CONFIG["references"] if section.get("references") else ""
        )
        description = description or section["description"].format_map(ctx)
        helptext = section["helptext"].format_map(ctx) + references
        content |= CONFIG["parent"] | {
            "section_name": section["section_name"],
            "description": description,
            "helptext": helptext,
        }
    pconfig = extra.pop(
        "pconfig", section.get("pconfig") if plot_type != "html" else None
    )
    content |= (
        {"plot_type": plot_type}
        | ({"pconfig": pconfig} if pconfig else {})
        | extra
        | {"data": data}
    )
    with open(f"{section_id}_mqc.json", "w") as handle:
        json.dump(
            content,
            handle,
            indent=4,
            allow_nan=False,
            default=lambda value: value.item(),
        )


def splices_text(count):
    return f"{count / 1e6:.2f} M" if count >= 1e5 else f"{count:,.0f}"


def contrasts_by_pair(contrasts):
    """The contrasts of each (treatment, control) pair, in contrastsheet order; contrasts can differ in blocking only."""
    by_pair = defaultdict(list)
    for name, pair in contrasts.items():
        by_pair[pair].append(name)
    return by_pair


class Reverses:
    def __init__(self, contrasts):
        self.order, by_pair = (
            {name: index for index, name in enumerate(contrasts)},
            contrasts_by_pair(contrasts),
        )
        self.other = {
            name: by_pair[pair[::-1]][0]
            for name, pair in contrasts.items()
            if pair[::-1] in by_pair
        }

    def is_reverse(self, contrast):
        return (
            contrast in self.other
            and self.order[contrast] > self.order[self.other[contrast]]
        )

    def label(self, contrast):
        return (
            f"{contrast} (reverse of {self.other[contrast]})"
            if self.is_reverse(contrast)
            else contrast
        )


class Genes:
    def __init__(self):
        self.names, self.ids_by_name = {}, {}

    def load(self, path):
        table = read_table(path)
        names = dict(
            zip(
                table["gene_id"].map(self.strip_version),
                table.get("gene_name", table["gene_id"]),
            )
        )
        self.names = {gene: name or gene for gene, name in names.items()}
        ids = defaultdict(set)
        for gene, name in self.names.items():
            ids[name].add(gene)
        self.ids_by_name = {
            name: genes.pop() for name, genes in ids.items() if len(genes) == 1
        }

    def strip_version(self, gene):
        stem, dot, suffix = gene.rpartition(".")
        return (
            stem
            if dot
            and suffix.isdigit()
            and (stem.startswith("ENS") or stem in self.names)
            else gene
        )

    def key(self, gene):
        gene = self.strip_version(str(gene).strip().strip('"'))
        return (
            None
            if gene in ("", "NA", ".")
            else gene
            if gene in self.names
            else self.ids_by_name.get(gene, gene)
        )

    def name(self, gene):
        return self.names.get(gene, gene)


#
# Parsers: each reads one input, records per-sample values on the run and returns feature frames
#

PARSERS = {}


def parses(*kinds):
    def register(function):
        PARSERS.update(dict.fromkeys(kinds, function))
        return function

    return register


@parses("gene_annotation")
def gene_annotation(run, label, path):
    run.genes.load(path)


@parses("samplesheet")
def samplesheet(run, label, path):
    run.condition_of |= dict(
        read_table(path, sep=",")[["sample", "condition"]].itertuples(index=False)
    )


@parses("contrasts")
def contrastsheet(run, label, path):
    table = read_table(path, sep=",")
    run.contrasts |= {
        row.contrast: Contrast(row.treatment, row.control) for row in table.itertuples()
    }
    present = set(run.condition_of.values())
    missing = [
        (name, condition)
        for name, pair in run.contrasts.items()
        for condition in pair
        if present and condition not in present
    ]
    if missing:
        raise ValueError(
            f"condition '{missing[0][1]}' of contrast '{missing[0][0]}' has no sample in the samplesheet"
        )


@parses("star_log")
def star_log(run, label, path):
    with open(path) as handle:
        values = {
            name.strip(): value.strip()
            for name, _, value in (line.partition("|") for line in handle)
        }
    total, annotated = numbers(
        pd.Series(
            [
                values.get(f"Number of splices: {name}")
                for name in ("Total", "Annotated (sjdb)")
            ]
        )
    ).fillna(0)
    if total:
        run.splices[label["id"]] = total
        run.general_stats[label["id"]] |= {
            "star_annotated_splices": 100 * annotated / total,
            "star_splices": total / 1e6,
        }


@parses("dexseq_counts")
def dexseq_counts(run, label, path):
    counts = pd.read_csv(
        path, sep=TAB, header=None, index_col=0, keep_default_na=False
    ).iloc[:, 0]
    unassigned = counts.index.isin(CONFIG["dexseq_counters"])
    in_bins = counts[
        ~unassigned & (counts.index != "") & ~counts.index.str.startswith("_")
    ].sum()
    total = in_bins + counts[unassigned].sum()
    run.general_stats[label["id"]]["dexseq_reads_in_bins"] = (
        100 * in_bins / total if total else None
    )


@parses("leafcutter_numers")
def leafcutter_numers(run, label, path):
    counts = pd.read_csv(path, sep=" ")
    for sample, reads in counts.items():
        detected = 100 * (reads > 0).sum() / len(reads) if len(reads) else None
        run.general_stats[sample] |= {
            "leafcutter_junction_reads": reads.sum(),
            "leafcutter_introns_detected": detected,
        }


@parses("leafcutter_strands")
def leafcutter_strands(run, label, path):
    table = read_table(path)
    for sample, value in zip(table["Sample"], numbers(table["Reads strandless (%)"])):
        run.general_stats[sample]["leafcutter_reads_strandless"] = (
            None if pd.isna(value) else value
        )


@parses("suppa_psi")
def suppa_psi(run, label, path):
    psi = (
        pd.read_csv(path, sep=TAB, index_col=0).apply(numbers)
        if os.path.getsize(path)
        else pd.DataFrame()
    )
    if len(psi) and any(psi.columns):
        # SUPPA2 releases before 2.0 wrote -1 for a missing PSI
        run.psi[label["route"]] = psi.where(psi >= 0)


def scan(run, label, path):
    spec = CONFIG["kinds"][label["kind"]]
    wanted, prefix = (
        {spec[key] for key in ("q", "p", "gene") if key in spec},
        spec.get("effect", ""),
    )
    table = read_table(
        path,
        sep=spec.get("sep", TAB),
        usecols=lambda column: (
            column in wanted or bool(prefix) and column.startswith(prefix)
        ),
    )
    qvalue = numbers(table[spec["q"]])
    table, qvalue = table[qvalue.notna()], qvalue.dropna()
    effect = (
        prefix
        if prefix in table
        else next(
            (column for column in table if prefix and column.startswith(prefix)), None
        )
    )
    columns = {
        "gene": table[spec["gene"]].map(run.genes.key) if "gene" in spec else None,
        "p": numbers(table[spec["p"]]) if "p" in spec else None,
    }
    significant = qvalue < run.args[spec.get("cutoff", "fdr")]
    return [
        features(
            spec["method"],
            spec["unit"],
            significant,
            effect=numbers(table[effect]) if effect else None,
            **columns,
        )
    ]


PARSERS |= dict.fromkeys(CONFIG["kinds"], scan)


@parses("stager")
def stager(run, label, path):
    """Genes from the screening stage, and transcripts from the confirmation stage, which only tests the transcripts of the genes that pass screening."""
    table, alpha = read_table(path), run.args["stager_alpha"]
    gene_padj = numbers(table["gene"])
    table, screened = table[gene_padj.notna()], gene_padj.dropna() < alpha
    passed = screened.groupby(table["geneID"], sort=False).last()
    confirmed = numbers(table.loc[screened, "transcript"]) < alpha
    method = "DEXSeq DTU + stageR"
    return [
        features(method, "genes", passed, gene=passed.index.map(run.genes.key)),
        features(method, "transcripts", confirmed),
    ]


@parses("suppa_events", "suppa_isoforms")
def suppa_diffsplice(run, label, path):
    table = read_table(path)
    dpsi, pvalue = numbers(table.iloc[:, 0]), numbers(table.iloc[:, 1])
    tested = dpsi.notna() & pvalue.notna()
    dpsi, pvalue, ids = dpsi[tested], pvalue[tested], dpsi.index[tested].to_series()
    events = label["kind"] == "suppa_events"
    types = {code: name for code, (name, *_) in CONFIG["suppa_event_types"].items()}
    columns = {
        "gene": ids.str.split(";").str[0].map(run.genes.key),
        "p": pvalue,
        "effect": dpsi,
    }
    if events:
        columns["event_type"] = (
            ids.str.split(";").str[1].str.split(":").str[0].map(types)
        )
    significant = (pvalue < run.args["suppa_alpha"]) & (
        dpsi.abs() > run.args["suppa_dpsi"]
    )
    return [
        features(
            "SUPPA2", "events" if events else "transcripts", significant, **columns
        )
    ]


@parses("isa_features")
def isoformswitchanalyzer(run, label, path):
    """Isoforms and genes with a switch q-value, for every contrast of the conditions of each comparison, as IsoformSwitchAnalyzeR ignores blocking."""
    table = read_table(path, sep=",")
    qvalue = numbers(table["isoform_switch_q_value"])
    table, qvalue = table[qvalue.notna()], qvalue.dropna()
    delta = numbers(table["dIF"])
    significant = (qvalue < run.args["isa_alpha"]) & (delta.abs() > run.args["isa_dif"])
    by_pair, method, frames = (
        contrasts_by_pair(run.contrasts),
        "IsoformSwitchAnalyzeR",
        [],
    )
    for pair, rows in table.groupby(
        ["condition_2", "condition_1"], sort=False
    ).groups.items():
        genes = table.loc[rows, "gene_id"].map(run.genes.key)
        switched = significant[rows].groupby(genes, sort=False, dropna=False).any()
        for contrast in by_pair.get(pair, ["-".join(pair)]):
            frames.append(
                features(
                    method,
                    "isoforms",
                    significant[rows],
                    effect=delta[rows],
                    contrast=contrast,
                )
            )
            frames.append(
                features(
                    method,
                    "genes with switches",
                    switched,
                    gene=switched.index,
                    contrast=contrast,
                )
            )
    return frames


@parses("leafcutter_ds")
def leafcutter_ds(run, label, path):
    table = read_table(path)
    padj, success = numbers(table["p.adjust"]), table["status"] == "Success"
    tested = success & padj.notna()
    significant = padj[tested] < run.args["fdr"]
    outcome = table["status"].map(CONFIG["leafcutter_skip_statuses"]).fillna("error")
    outcome[success] = "no_pvalue"
    outcome[tested] = significant.map({True: "significant", False: "not_significant"})
    run.leafcutter[label["id"]] = dict(Counter(outcome))
    genes = (
        [
            [gene for gene in map(run.genes.key, names.split(",")) if gene]
            for names in table.loc[tested, "genes"]
        ]
        if "genes" in table
        else None
    )
    return [
        features(
            "LeafCutter",
            "clusters",
            significant,
            gene=genes,
            p=numbers(table.loc[tested, "p"]) if "p" in table else None,
        )
    ]


@parses("rmats_summary")
def rmats_summary(run, label, path):
    table = read_table(path)
    tested, significant = (
        table[column].astype(int).sum()
        for column in ("TotalEventsJC", "SignificantEventsJC")
    )
    row = ("rMATS", label["route"], "JC events", label["id"], tested, significant)
    run.counts.append(dict(zip([*GROUP, "tested", "significant"], row)))


@parses("rmats_jc")
def rmats_jc(run, label, path):
    """Genes, p-values and effect sizes of the rMATS events; their numbers come from the rMATS summary."""
    table = read_table(
        path,
        usecols=lambda column: (
            column in ("GeneID", "PValue", "FDR", "IncLevelDifference")
        ),
    )
    fdr = numbers(table["FDR"])
    table, significant = table[fdr.notna()], fdr.dropna() <= LIMITS["rmats_fdr"]
    columns = {
        "gene": table["GeneID"].map(run.genes.key),
        "p": numbers(table["PValue"]),
        "effect": numbers(table["IncLevelDifference"]),
    }
    return [features("rMATS", "JC events", significant, counted=False, **columns)]


def unknown_kind(run, label, path):
    raise ValueError(f"Unknown splicing result kind '{label['kind']}' for {path}")


@contextmanager
def fail_soft(run, kind, path):
    try:
        yield
    except Exception as error:  # noqa: BLE001 any input error leaves that input out, not the report
        failure = Failure(
            CONFIG["inputs"].get(kind, {}).get("name", kind),
            os.path.basename(path),
            f"{type(error).__name__}: {error}",
        )
        print(
            f"WARNING: could not summarise {failure.name} ({failure.filename}): {failure.error}. Its results are left out of the report.",
            file=sys.stderr,
        )
        run.failures.append(failure)


#
# Sections
#


def write_general_stats(run):
    headers = CONFIG["general_stats"]
    columns = [
        column
        for column in headers
        if any(values.get(column) is not None for values in run.general_stats.values())
    ]
    data = {
        sample: {
            column: values[column]
            for column in columns
            if values.get(column) is not None
        }
        for sample, values in sorted(run.general_stats.items())
    }
    emit(
        "splicing_sample_metrics",
        "generalstats",
        data,
        headers={column: headers[column] for column in columns},
    )


def leafcutter_settings(args):
    defaults = CONFIG["leafcutter_defaults"]
    settings = {
        key: args[f"leafcutter_min_{name}"]
        for key, name in (
            ("g", "samples_per_group"),
            ("i", "samples_per_intron"),
            ("c", "coverage"),
        )
    }
    flags = "-g {g} -i {i} -c {c}"
    return (
        settings,
        flags.format(**settings),
        flags.format(**defaults),
        any(settings[key] < defaults[key] for key in defaults),
    )


def write_verdict(run):
    section, sizes = "splicing_verdict", run.sizes()
    say = lambda key, **values: message(section, key, **values)
    if not sizes and not run.splices and not run.failures:
        return
    facts, evidence = [], []
    if sizes:
        designs = {f"{treatment} vs {control}" for treatment, control in sizes.values()}
        listed = ", ".join(
            f"{treatment} vs {control} ({name})"
            for name, (treatment, control) in sorted(sizes.items())
        )
        facts.append(
            say("design", designs=designs.pop())
            if len(designs) == 1
            else say("designs", designs=listed)
        )
    _, flags, defaults, below = leafcutter_settings(run.args)
    if "LeafCutter" in run.ran and below:
        facts.append(say("leafcutter", flags=flags, defaults=defaults))
    if run.splices:
        values = pd.Series(run.splices)
        spread = splices_text(values.min()) + (
            "" if values.min() == values.max() else f" to {splices_text(values.max())}"
        )
        facts.append(say("splices", median=splices_text(values.median()), range=spread))
    smallest = min((size for pair in sizes.values() for size in pair), default=None)
    if smallest is not None and run.ran & DTU and smallest <= 2:
        evidence.append(say("evidence_dtu"))
    if smallest is not None and "LeafCutter" in run.ran and smallest < 4:
        evidence.append(say("evidence_leafcutter"))
    text = " ".join(facts)
    if evidence:
        sentence = "; ".join(evidence)
        text += f" {sentence[0].upper()}{sentence[1:]}. " + say("exploratory")
    content = say("facts", message=text) if text else ""
    content += "".join(
        say(
            "failure",
            **{key: html.escape(value) for key, value in failure._asdict().items()},
        )
        for failure in run.failures
    )
    emit(section, "html", content)


def design_problems(run, n_treatment, n_control):
    say = lambda key, **values: message(
        "splicing_design_check",
        key,
        n_treatment=n_treatment,
        n_control=n_control,
        **values,
    )
    args, ran, smallest, pair = (
        run.args,
        run.ran,
        min(n_treatment, n_control),
        (n_treatment, n_control),
    )
    cannot, below, evidence = [], [], []
    variance = [name for analysis, name in VARIANCE.items() if analysis in ran]
    if smallest < 2 and variance:
        cannot.append(say("variance", names=", ".join(variance)))
    for analysis, option, key in (
        ("SUPPA2", "suppa_paired", "suppa_paired"),
        ("rMATS", "rmats_paired_stats", "rmats_paired"),
    ):
        if analysis in ran and args[option] and n_treatment != n_control:
            cannot.append(say(key))
    if "LeafCutter" in ran:
        settings, flags, defaults, below_defaults = leafcutter_settings(args)
        minimum = args["leafcutter_min_samples_per_condition"]
        cannot += (
            [say("leafcutter_minimum", minimum=minimum)] if smallest < minimum else []
        )
        below += (
            [say("leafcutter_defaults", flags=flags, defaults=defaults)]
            if below_defaults
            else []
        )
        if smallest < 4:
            below.append(say("leafcutter_few"))
            evidence.append(say("evidence_leafcutter_few"))
        if pair == (2, 2) and settings["i"] < 4:
            below.append(say("leafcutter_i", i=settings["i"]))
            evidence.append(say("evidence_leafcutter_i", i=settings["i"]))
    needed, total = (
        max(
            args[f"min_samps_{name}"]
            for name in ("gene_expr", "feature_expr", "feature_prop")
        ),
        len(run.condition_of),
    )
    if "DRIMSeq" in ran and total < needed:
        cannot.append(say("drimseq", needed=needed, total=total))
    if ran & DTU and smallest <= 3:
        evidence.append(say("evidence_dtu_2" if smallest <= 2 else "evidence_dtu_3"))
    if pair == (2, 2):
        evidence.append(say("evidence_permutation"))
    return cannot, below, evidence


def write_design_check(run):
    section = "splicing_design_check"
    if not run.contrasts or not run.condition_of:
        return
    leafcutter, drimseq = "LeafCutter" in run.ran, "DRIMSeq" in run.ran
    shown = {"splices": bool(run.splices), "drimseq_gene": drimseq} | dict.fromkeys(
        ("leafcutter_g", "leafcutter_i", "leafcutter_c"), leafcutter
    )
    data = {}
    for name, (n_treatment, n_control) in sorted(run.sizes().items()):
        contrast = run.contrasts[name]
        row = {
            "treatment": contrast.treatment,
            "control": contrast.control,
            "n_treatment": n_treatment,
            "n_control": n_control,
        }
        row |= (
            {
                f"leafcutter_{key}": value
                for key, value in leafcutter_settings(run.args)[0].items()
            }
            if leafcutter
            else {}
        )
        row |= {"drimseq_gene": run.args["min_samps_gene_expr"]} if drimseq else {}
        splices = [
            run.splices[sample]
            for sample in run.samples(contrast.treatment)
            | run.samples(contrast.control)
            if sample in run.splices
        ]
        row |= {"splices": np.median(splices) / 1e6} if splices else {}
        cannot, below, evidence = design_problems(run, n_treatment, n_control)
        row["status"] = message(
            section, "cannot" if cannot else "below" if below else "meets"
        )
        details = [
            message(section, key, items="; ".join(items))
            for key, items in (("details_cannot", cannot), ("details_below", below))
            if items
        ]
        data[name] = row | {
            "details": " ".join(details),
            "evidence": "; ".join(evidence),
        }
    emit(
        section,
        "table",
        data,
        headers={
            key: header
            for key, header in SECTIONS[section]["headers"].items()
            if shown.get(key, True)
        },
    )


def write_concordance(run):
    section, route = "suppa_psi_concordance", first_route(run.psi)
    if route is None:
        return
    psi = run.psi[route]
    complete = psi.dropna()
    informative = complete[complete.mean(axis=1).between(0.05, 0.95)]
    if psi.shape[1] < 2 or len(informative) < LIMITS["psi_concordance_min_events"]:
        return
    samples = sorted(
        psi.columns, key=lambda sample: (run.condition_of.get(sample, ""), sample)
    )
    correlation = informative[samples].corr(method="spearman").round(3).astype(object)
    names = [
        f"{sample} ({run.condition_of[sample]})"
        if sample in run.condition_of
        else sample
        for sample in samples
    ]
    few_events = LIMITS["psi_concordance_few_events"]
    description = SECTIONS[section]["description"].format(
        total=f"{len(informative):,}",
        route=with_route("", route, SECTIONS[section]["messages"]["route"])
        if len(run.psi) > 1
        else "",
        few=message(section, "few", few_events=few_events)
        if len(informative) < few_events
        else "",
    )
    matrix = correlation.where(correlation.notna(), None).to_numpy().tolist()
    emit(
        section,
        "heatmap",
        matrix,
        {"few_events": few_events},
        description,
        xcats=names,
        ycats=names,
    )


def overview(run, calls):
    counts = (
        calls[calls["counted"]]
        .groupby(GROUP, sort=False)["significant"]
        .agg(tested="size", significant="sum")
        .reset_index()
    )
    counts = pd.concat(
        [counts, pd.DataFrame(run.counts, columns=[*GROUP, "tested", "significant"])],
        ignore_index=True,
    )
    significant = calls[calls["significant"] & calls["effect"].notna()]
    effects = (
        significant["effect"]
        .astype(float)
        .abs()
        .groupby([significant[column] for column in GROUP])
        .median()
    )
    sources = [
        (
            row.method,
            row.route,
            UNITS[row.method][row.unit].get("effect_from", row.unit),
            row.contrast,
        )
        for row in counts.itertuples()
    ]
    counts["effect"] = [
        effects.get(source) if "effect_unit" in UNITS[row.method][row.unit] else None
        for source, row in zip(sources, counts.itertuples())
    ]
    rank = {"method": list(METHODS).index, "route": ROUTES.index}
    return counts.sort_values(
        GROUP,
        key=lambda column: (
            column.map(rank[column.name]) if column.name in rank else column
        ),
        kind="stable",
    )


def column_headers(run, row, templates, messages):
    unit = UNITS[row.method][row.unit]
    ctx = {
        "method": row.method,
        "unit": row.unit,
        "title": f"{row.method} {row.unit}"
        + with_route("", row.route, messages["route"]),
        "route_text": with_route("", row.route, messages["route_text"]),
        "threshold": unit["threshold"].format_map(run.args),
        "effect": unit.get("effect", CONFIG["no_effect_cutoff"]).format_map(run.args),
        "tested_text": unit.get("tested") or messages["tested"].format(unit=row.unit),
        "effect_unit": unit.get("effect_unit"),
    }
    return {
        name: {
            field: value.format_map(ctx) if field in ("title", "description") else value
            for field, value in templates[name].items()
        }
        for name in templates
    }


def write_differential(run, counts, reverses):
    section = "differential_splicing"
    templates, messages = SECTIONS[section]["headers"], SECTIONS[section]["messages"]
    primary_routes = counts.groupby("method")["route"].agg(
        lambda routes: first_route(set(routes))
    )
    contrasts = sorted(set(counts["contrast"]))
    headers, data = {}, defaultdict(dict)
    if contrasts_reversed := set(contrasts) & reverses.other.keys():
        headers["reciprocal"] = templates["reciprocal"]
        for contrast in contrasts_reversed:
            data[contrast]["reciprocal"] = reverses.other[contrast]
    for row in counts.itertuples():
        primary = row.route == primary_routes[row.method]
        key = "_".join(
            "".join(
                char if char.isalnum() else " "
                for char in f"{row.method} {row.route} {row.unit}".lower()
            ).split()
        )
        if key not in headers:
            column = column_headers(run, row, templates, messages)
            headers[key] = column["significant"] | {
                "hidden": not (
                    primary and METHODS[row.method].get("primary") == row.unit
                )
            }
            headers[f"{key}_tested"] = column["tested"]
            if "effect_unit" in UNITS[row.method][row.unit] and primary:
                headers[f"{key}_effect"] = column["effect"]
        data[row.contrast] |= {key: row.significant, f"{key}_tested": row.tested}
        data[row.contrast] |= (
            {f"{key}_effect": row.effect} if primary and pd.notna(row.effect) else {}
        )
    ctx = {"reciprocal": messages["reciprocal"] if "reciprocal" in headers else ""}
    emit(
        section,
        "table",
        {contrast: data[contrast] for contrast in contrasts},
        ctx,
        headers=headers,
    )


def write_leafcutter(run, reverses):
    section = "leafcutter_cluster_test_outcomes"
    categories = {
        key: category | {"name": category["name"].format_map(run.args)}
        for key, category in SECTIONS[section]["categories"].items()
    }
    data = {
        reverses.label(contrast): counts
        for contrast, counts in sorted(run.leafcutter.items())
    }
    emit(section, "bargraph", data, run.args, categories=categories)


def write_suppa_events(run, calls, reverses):
    section, minimum = "suppa_significant_events", LIMITS["suppa_type_min_significant"]
    events = calls[(calls["method"] == "SUPPA2") & (calls["unit"] == "events")]
    route = first_route(set(events["route"]))
    if route is None:
        return
    several = events["route"].nunique() > 1
    ctx = {
        "min_significant": minimum,
        "significance": "{suppa_correction} p < {suppa_alpha:g}, |dPSI| > {suppa_dpsi:g}".format_map(
            run.args
        ),
        "route": with_route("", route, SECTIONS[section]["messages"]["route"])
        if several
        else "",
    }
    events = events[events["route"] == route]
    found = events.groupby("contrast")["significant"].any()
    description = SECTIONS[section]["description"].format_map(ctx)
    if not found.any():
        return emit(section, "html", message(section, "none"), ctx, description)
    if not found.all():
        description += message(
            section,
            "without",
            contrasts=", ".join(map(reverses.label, found.index[~found])),
        )
    typed = events[events["event_type"].notna()]
    category = lambda event_type, higher: message(
        section,
        "category",
        event_type=event_type,
        direction="higher" if higher else "lower",
    )
    counts, percentages = {}, {}
    for contrast in found.index[found]:
        tested = typed[typed["contrast"] == contrast]
        significant = tested[tested["significant"]]
        pairs = Counter(
            zip(significant["event_type"], significant["effect"].astype(float) > 0)
        )
        per_type, of_type = (
            significant["event_type"].value_counts(),
            tested["event_type"].value_counts(),
        )
        counts[reverses.label(contrast)] = {
            category(*pair): count for pair, count in pairs.items()
        }
        shares = {
            category(*pair): 100 * count / of_type[pair[0]]
            for pair, count in pairs.items()
            if per_type[pair[0]] >= minimum
        }
        percentages |= {reverses.label(contrast): shares} if shares else {}
    colours = CONFIG["suppa_event_types"].values()
    categories = {
        category(name, higher): {"name": category(name, higher), "color": colour}
        for name, *shades in colours
        for higher, colour in zip((True, False), shades)
    }
    pconfig = dict(SECTIONS[section]["pconfig"])
    if percentages:
        pconfig["data_labels"] = [
            {key: value.format_map(ctx) for key, value in label.items()}
            for label in SECTIONS[section]["data_labels"]
        ]
        return emit(
            section,
            "bargraph",
            [counts, percentages],
            ctx,
            description,
            pconfig=pconfig,
            categories=[categories] * 2,
        )
    emit(
        section,
        "bargraph",
        counts,
        ctx,
        description,
        pconfig=pconfig,
        categories=categories,
    )


def write_agreement(run, calls, reverses):
    section, maximum = "splicing_method_agreement", LIMITS["agreement_max_genes"]
    primary = calls[
        calls["unit"]
        == calls["method"].map(lambda method: METHODS[method].get("primary"))
    ]
    calling = primary[~primary["contrast"].map(reverses.is_reverse).astype(bool)]
    calling = calling[calling.groupby("contrast")["method"].transform("nunique") >= 2]
    if calling.empty:
        return
    two_routes = calling.groupby(["contrast", "method"])["route"].nunique().gt(1).any()
    hits = calling[calling["significant"]].explode("gene").dropna(subset=["gene"])
    routes = hits.groupby(["contrast", "gene", "method"])["route"].nunique()
    rows, capped = defaultdict(list), []
    for (contrast, gene), called in routes.groupby(level=["contrast", "gene"]):
        found = called.droplevel(["contrast", "gene"])
        names = [method for method in METHODS if method in found.index]
        families = sorted({METHODS[method]["family"] for method in names})
        row = {
            "gene_id": gene,
            "n_families": len(families),
            "n_methods": len(names),
            "families": ", ".join(families),
            "methods": ", ".join(names),
        }
        row |= (
            {
                "both_routes": ", ".join(
                    method for method in names if found[method] > 1
                )
                or "No"
            }
            if two_routes
            else {}
        )
        rows[contrast] += [row] if len(families) >= 2 else []
    table = {}
    for contrast, genes in sorted(rows.items()):
        genes.sort(
            key=lambda row: (
                -row["n_families"],
                -row["n_methods"],
                run.genes.name(row["gene_id"]),
            )
        )
        capped += [f"{contrast} ({len(genes):,})"] if len(genes) > maximum else []
        table |= {
            f"{run.genes.name(row['gene_id'])} ({contrast})": row
            for row in genes[:maximum]
        }
    description = SECTIONS[section]["description"]
    reverse_items = [
        message(section, "reverse_item", contrast=contrast, other=other)
        for contrast, other in sorted(reverses.other.items())
        if reverses.is_reverse(contrast) and contrast in set(primary["contrast"])
    ]
    description += (
        message(section, "reverse", contrasts=", ".join(reverse_items))
        if reverse_items
        else ""
    )
    description += (
        message(section, "capped", maximum=maximum, contrasts=", ".join(capped))
        if capped
        else ""
    )
    if not table:
        return emit(section, "html", message(section, "none"), description=description)
    headers = {
        key: header
        for key, header in SECTIONS[section]["headers"].items()
        if key != "both_routes" or two_routes
    }
    emit(section, "table", table, description=description, headers=headers)


def write_pvalues(calls, reverses):
    section, bins = "splicing_pvalue_histograms", LIMITS["p_bins"]
    keys = list(zip(calls["method"], calls["unit"], calls["route"]))
    labels = {
        key: with_route(UNITS[key[0]][key[1]]["histogram"], key[2])
        for key in set(keys)
        if "histogram" in UNITS[key[0]][key[1]]
    }
    pvalues = calls.assign(
        series=[labels.get(key) for key in keys], p=numbers(calls["p"])
    )
    pvalues = pvalues[pvalues["series"].notna() & pvalues["p"].between(0, 1)]
    pvalues = pvalues[
        pvalues.groupby(["contrast", "series"], sort=False)["p"].transform("size")
        >= LIMITS["p_min_features"]
    ]
    contrasts = sorted(set(pvalues["contrast"]))
    if not contrasts:
        return

    def histogram(p):
        counts = np.bincount(
            np.minimum((p * bins).astype(int), bins - 1), minlength=bins
        )
        return {
            f"{(position + 0.5) / bins:.3f}": 100 * count / len(p)
            for position, count in enumerate(counts)
        }

    data = [
        {
            series: histogram(p)
            for series, p in pvalues[pvalues["contrast"] == contrast].groupby(
                "series", sort=False
            )["p"]
        }
        for contrast in contrasts
    ]
    pconfig = SECTIONS[section]["pconfig"] | {
        "data_labels": [
            {"name": reverses.label(contrast), "ylab": "% of tested features"}
            for contrast in contrasts
        ]
    }
    emit(
        section,
        "linegraph",
        data,
        {
            "minimum": LIMITS["p_min_features"],
            "bins": bins,
            "per_bin": f"{100 / bins:g}",
        },
        pconfig=pconfig,
    )


def main():
    run = Run(parse_args(EXT_ARGS))
    for label, path in sorted(
        zip(LABELS, STAGED),
        key=lambda item: (
            item[0]["kind"] not in ("gene_annotation", "samplesheet", "contrasts")
        ),
    ):
        with fail_soft(run, label["kind"], path):
            frames = PARSERS.get(label["kind"], unknown_kind)(run, label, path) or []
            run.frames += [
                frame.assign(
                    route=label["route"], contrast=frame.get("contrast", label["id"])
                )
                for frame in frames
            ]
            run.ran.add(CONFIG["inputs"].get(label["kind"], {}).get("analysis"))
    psi = run.psi.get(first_route(run.psi))
    for sample in psi if psi is not None else []:
        run.general_stats[sample]["suppa_events_quantified"] = (
            100 * psi[sample].notna().sum() / len(psi)
        )

    reverses = Reverses(run.contrasts)
    columns = [*GROUP, "significant", "gene", "p", "effect", "event_type", "counted"]
    frames = [frame.reindex(columns=columns) for frame in run.frames] or [
        pd.DataFrame(columns=columns)
    ]
    calls = pd.concat(frames, ignore_index=True).astype(
        {"counted": bool, "significant": bool}
    )
    if run.general_stats:
        write_general_stats(run)
    write_verdict(run)
    write_design_check(run)
    write_concordance(run)
    if len(counts := overview(run, calls)):
        write_differential(run, counts, reverses)
    if run.leafcutter:
        write_leafcutter(run, reverses)
    write_suppa_events(run, calls, reverses)
    write_agreement(run, calls, reverses)
    write_pvalues(calls, reverses)
    versions = {"python": platform.python_version()} | {
        package.__name__: package.__version__ for package in (np, pd, yaml)
    }
    with open("versions.yml", "w") as handle:
        yaml.safe_dump({PROCESS: versions}, handle)


if __name__ == "__main__":
    main()
