#!/usr/bin/env python

# Summarise the result tables of the splicing analyses into MultiQC custom-content files:
# General Statistics columns, and the sections of a Splicing parent section, whose text,
# headers and tool tables come from the section config.

import argparse
import html
import json
import os
import platform
import re
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
SECTIONS, LIMITS, METHODS, UNITS, KINDS = (
    CONFIG[key] for key in ("sections", "limits", "methods", "units", "kinds")
)
ROUTES = ["STAR", "Salmon", ""]
DTU = {"DRIMSeq", "IsoformSwitchAnalyzeR"}
GROUP = ["method", "route", "unit", "contrast"]
# Read first, as the contrastsheet is checked against the samplesheet and the others use both
PRELOAD = ["gene_annotation", "samplesheet", "contrasts"]
EVENT_TYPES = {code: name for code, (name, *_) in CONFIG["suppa_event_types"].items()}
Contrast = namedtuple("Contrast", "treatment control")
Failure = namedtuple("Failure", "name filename error")


def parse_args():
    parser = argparse.ArgumentParser()
    for name, default in CONFIG["args"].items():
        flag = isinstance(default, bool)
        parser.add_argument(
            f"--{name}",
            type=(lambda value: value.lower() == "true") if flag else type(default),
            default=default,
        )
    args = vars(parser.parse_args(shlex.split(EXT_ARGS))) | LIMITS
    leafcutter = {flag: args[name] for flag, name in CONFIG["leafcutter_flags"].items()}
    return args | {
        "suppa_correction": "gene-corrected"
        if args["suppa_gene_correction"]
        else "uncorrected",
        "flags": "-g {g} -i {i} -c {c}".format(**leafcutter),
        "leafcutter_default_flags": CONFIG["leafcutter_default_flags"],
        "leafcutter_below_defaults": any(
            leafcutter[flag] < default
            for flag, default in CONFIG["leafcutter_defaults"].items()
        ),
    }


ARGS = parse_args()


class Genes:
    """Gene IDs without version suffixes, with the names some tools report mapped to their ID."""

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
        counts = Counter(self.names.values())
        self.ids_by_name = {
            name: gene for gene, name in self.names.items() if counts[name] == 1
        }

    def strip_version(self, gene):
        stem, dot, suffix = gene.rpartition(".")
        if dot and suffix.isdigit() and (stem.startswith("ENS") or stem in self.names):
            return stem
        return gene

    def key(self, gene):
        gene = self.strip_version(str(gene).strip().strip('"'))
        if gene in ("", "NA", "."):
            return None
        return gene if gene in self.names else self.ids_by_name.get(gene, gene)

    def name(self, gene):
        return self.names.get(gene, gene)


@dataclass
class Run:
    genes: Genes = field(default_factory=Genes)
    contrasts: dict = field(default_factory=dict)
    condition_of: dict = field(default_factory=dict)
    frames: list = field(default_factory=list)
    general_stats: dict = field(default_factory=lambda: defaultdict(dict))
    splices: dict = field(default_factory=dict)
    psi: dict = field(default_factory=dict)
    leafcutter: dict = field(default_factory=dict)
    failures: list = field(default_factory=list)
    ran: set = field(default_factory=set)

    def samples(self, condition):
        return {sample for sample, of in self.condition_of.items() if of == condition}

    def sizes(self):
        """Samples of the treatment and control of each contrast whose conditions are in the samplesheet."""
        counts = Counter(self.condition_of.values())
        return {
            name: (counts[treatment], counts[control])
            for name, (treatment, control) in self.contrasts.items()
            if treatment in counts and control in counts
        }


class Reverses:
    """Contrasts that test the same comparison as another; the later of each pair in the contrastsheet is the reverse of the first."""

    def __init__(self, contrasts):
        names, by_pair = list(contrasts), contrasts_by_pair(contrasts)
        self.other = {
            name: by_pair[pair[::-1]][0]
            for name, pair in contrasts.items()
            if pair[::-1] in by_pair
        }
        self.later = {
            name
            for name, other in self.other.items()
            if names.index(name) > names.index(other)
        }

    def label(self, contrast):
        return (
            f"{contrast} (reverse of {self.other[contrast]})"
            if contrast in self.later
            else contrast
        )


def read_table(path, sep=TAB, **kwargs):
    """A table as strings, as written, so that numbers() treats every missing value alike."""
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, **kwargs)


def numbers(series):
    return pd.to_numeric(series, errors="coerce")


def tested(table, column):
    values = numbers(table[column])
    return table[values.notna()], values.dropna()


def features(method, unit, significant, counted=True, **columns):
    """Tested features of one analysis, one row each, with their gene, p-value and effect size where the table has them."""
    frame = pd.DataFrame({"significant": significant, **columns})
    return frame.assign(method=method, unit=unit, counted=counted).reset_index(
        drop=True
    )


def routed(message, route):
    return (
        CONFIG["route_messages"][message].format(route=CONFIG["routes"][route])
        if route
        else ""
    )


def first_route(routes):
    """The STAR and Salmon routes quantify the same samples, so per-sample summaries show one."""
    return next((route for route in ROUTES if route in routes), None)


def say(section, key, **values):
    return SECTIONS[section]["messages"][key].format_map(ARGS | values)


def dump(name, content):
    with open(f"{name}_mqc.json", "w") as handle:
        json.dump(
            content,
            handle,
            indent=4,
            allow_nan=False,
            default=lambda value: value.item(),
        )


def emit(section_id, data, ctx=None, plot_type=None, **extra):
    section, ctx = SECTIONS[section_id], ARGS | (ctx or {})
    plot_type = plot_type or section["plot_type"]
    helptext = section["helptext"].format_map(ctx)
    if section.get("references"):
        helptext += "\\n\\n" + CONFIG["references"]
    content = {"id": section_id} | CONFIG["parent"]
    content |= {
        "section_name": section["section_name"],
        "description": section["description"].format_map(ctx),
        "helptext": helptext,
        "plot_type": plot_type,
    }
    pconfig = extra.pop("pconfig", section.get("pconfig"))
    if pconfig and plot_type != "html":
        content["pconfig"] = pconfig
    dump(section_id, content | extra | {"data": data})


def splices_text(count):
    return f"{count / 1e6:.2f} M" if count >= 1e5 else f"{count:,.0f}"


def contrasts_by_pair(contrasts):
    """The contrasts of each (treatment, control) pair, in contrastsheet order; contrasts can differ in blocking only."""
    by_pair = defaultdict(list)
    for name, pair in contrasts.items():
        by_pair[pair].append(name)
    return by_pair


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
    contrasts = {
        row.contrast: Contrast(row.treatment, row.control)
        for row in read_table(path, sep=",").itertuples()
    }
    present = set(run.condition_of.values())
    for name, pair in contrasts.items():
        for condition in pair:
            if present and condition not in present:
                raise ValueError(
                    f"condition '{condition}' of contrast '{name}' has no sample in the samplesheet"
                )
    run.contrasts |= contrasts


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
    in_bins = counts[(counts.index != "") & ~counts.index.str.startswith("_")].sum()
    total = in_bins + counts[counts.index.isin(CONFIG["dexseq_counters"])].sum()
    if total:
        run.general_stats[label["id"]]["dexseq_reads_in_bins"] = 100 * in_bins / total


@parses("leafcutter_numers")
def leafcutter_numers(run, label, path):
    for sample, reads in pd.read_csv(path, sep=" ").items():
        run.general_stats[sample]["leafcutter_junction_reads"] = reads.sum()
        if len(reads):
            run.general_stats[sample]["leafcutter_introns_detected"] = (
                100 * (reads > 0).sum() / len(reads)
            )


@parses("leafcutter_strands")
def leafcutter_strands(run, label, path):
    table = read_table(path)
    for sample, value in zip(table["Sample"], numbers(table["Reads strandless (%)"])):
        run.general_stats[sample]["leafcutter_reads_strandless"] = value


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


def effect_column(table, prefix, contrast):
    """The prefix column, or the column of the contrast's treatment and control, as DEXSeq writes one per pair of levels."""
    if prefix in table:
        return prefix
    columns = [column for column in table if column.startswith(prefix)]
    if contrast and f"{prefix}{contrast.treatment}_{contrast.control}" in columns:
        return f"{prefix}{contrast.treatment}_{contrast.control}"
    return columns[0] if len(columns) == 1 else None


def scan(run, label, path):
    """A result table where a feature is tested when it has a q-value, as described by its kinds entry."""
    spec = KINDS[label["kind"]]
    wanted, prefix = (
        {spec[key] for key in ("q", "p", "gene") if key in spec},
        spec.get("effect"),
    )
    table = read_table(
        path,
        sep=spec.get("sep", TAB),
        usecols=lambda column: (
            column in wanted or bool(prefix) and column.startswith(prefix)
        ),
    )
    table, qvalue = tested(table, spec["q"])
    cutoff = ARGS[spec.get("cutoff", "fdr")]
    significant = qvalue <= cutoff if spec.get("inclusive") else qvalue < cutoff
    columns = {"gene": table[spec["gene"]]} if "gene" in spec else {}
    if "p" in spec:
        columns["p"] = numbers(table[spec["p"]])
    effect = (
        effect_column(table, prefix, run.contrasts.get(label["id"])) if prefix else None
    )
    if effect:
        columns["effect"] = numbers(table[effect])
    return [
        features(
            spec["method"],
            spec["unit"],
            significant,
            spec.get("counted", True),
            **columns,
        )
    ]


PARSERS |= {kind: scan for kind, spec in KINDS.items() if "q" in spec}


@parses("stager")
def stager(run, label, path):
    """Genes from the screening stage, and transcripts from the confirmation stage, which only tests the transcripts of the genes that pass screening."""
    table, gene_padj = tested(read_table(path), "gene")
    screened = gene_padj < ARGS["stager_alpha"]
    passed = screened.groupby(table["geneID"], sort=False).last()
    confirmed = numbers(table.loc[screened, "transcript"]) < ARGS["stager_alpha"]
    method = "DEXSeq DTU + stageR"
    return [
        features(method, "genes", passed, gene=passed.index),
        features(method, "transcripts", confirmed),
    ]


@parses("suppa_events", "suppa_isoforms")
def suppa_diffsplice(run, label, path):
    table = read_table(path)
    dpsi, pvalue = (
        numbers(table[next(c for c in table if c.endswith(suffix))])
        for suffix in ("_dPSI", "_p-val")
    )
    dpsi, pvalue = (
        dpsi[dpsi.notna() & pvalue.notna()],
        pvalue[dpsi.notna() & pvalue.notna()],
    )
    ids = dpsi.index.to_series().str.extract("^([^;]*);([^:]*)")
    columns = {"gene": ids[0], "p": pvalue, "effect": dpsi}
    if label["kind"] == "suppa_events":
        columns["event_type"] = ids[1].map(EVENT_TYPES)
    significant = (pvalue < ARGS["suppa_alpha"]) & (dpsi.abs() > ARGS["suppa_dpsi"])
    return [
        features(
            "SUPPA2",
            "events" if label["kind"] == "suppa_events" else "transcripts",
            significant,
            **columns,
        )
    ]


@parses("isa_features")
def isoformswitchanalyzer(run, label, path):
    """Isoforms and genes with a switch q-value, for every contrast of the conditions of each comparison, as IsoformSwitchAnalyzeR ignores blocking."""
    table, qvalue = tested(read_table(path, sep=","), "isoform_switch_q_value")
    delta = numbers(table["dIF"])
    significant = (qvalue < ARGS["isa_alpha"]) & (delta.abs() > ARGS["isa_dif"])
    by_pair, method, frames = (
        contrasts_by_pair(run.contrasts),
        "IsoformSwitchAnalyzeR",
        [],
    )
    for pair, rows in table.groupby(
        ["condition_2", "condition_1"], sort=False
    ).groups.items():
        switched = (
            significant[rows]
            .groupby(
                table.loc[rows, "gene_id"].map(run.genes.key), sort=False, dropna=False
            )
            .any()
        )
        isoforms = features(method, "isoforms", significant[rows], effect=delta[rows])
        switches = features(
            method, "genes with switches", switched, gene=switched.index
        )
        for contrast in by_pair.get(pair, ["-".join(pair)]):
            frames += [
                isoforms.assign(contrast=contrast),
                switches.assign(contrast=contrast),
            ]
    return frames


@parses("leafcutter_ds")
def leafcutter_ds(run, label, path):
    table = read_table(path)
    padj, success = numbers(table["p.adjust"]), table["status"] == "Success"
    keep = success & padj.notna()
    significant = padj[keep] < ARGS["fdr"]
    outcome = table["status"].map(CONFIG["leafcutter_skip_statuses"]).fillna("error")
    outcome[success] = "no_pvalue"
    outcome[keep] = significant.map({True: "significant", False: "not_significant"})
    run.leafcutter[label["id"]] = dict(Counter(outcome))
    columns = {"gene": table.loc[keep, "genes"]} if "genes" in table else {}
    if "p" in table:
        columns["p"] = numbers(table.loc[keep, "p"])
    return [features("LeafCutter", "clusters", significant, **columns)]


@parses("rmats_summary")
def rmats_summary(run, label, path):
    """The numbers of rMATS events, which the JC tables give the genes, p-values and effect sizes of."""
    table = read_table(path)
    total, significant = (
        table[column].astype(int).sum()
        for column in ("TotalEventsJC", "SignificantEventsJC")
    )
    return [
        features(
            "rMATS", "JC events", [True] * significant + [False] * (total - significant)
        )
    ]


@contextmanager
def fail_soft(run, name, filename):
    """Record an input or section that fails instead of stopping the report."""
    try:
        yield
    except Exception as error:  # noqa: BLE001 any failure leaves out its input or section, not the report
        failure = Failure(name, filename, f"{type(error).__name__}: {error}")
        print(
            f"WARNING: could not summarise {failure.name} ({failure.filename}): {failure.error}. Its results are left out of the report.",
            file=sys.stderr,
        )
        run.failures.append(failure)


#
# Sections: each writer takes its section id, the run, the tested features and their counts
#


def write_general_stats(run):
    data = {
        sample: {column: value for column, value in values.items() if pd.notna(value)}
        for sample, values in sorted(run.general_stats.items())
    }
    headers = {
        column: header
        for column, header in CONFIG["general_stats"].items()
        if any(column in values for values in data.values())
    }
    data = {
        sample: {column: values[column] for column in headers if column in values}
        for sample, values in data.items()
    }
    dump(
        "splicing_sample_metrics",
        {
            "id": "splicing_sample_metrics",
            "plot_type": "generalstats",
            "headers": headers,
            "data": data,
        },
    )


def write_verdict(section, run, calls, counts):
    sizes = run.sizes()
    if not sizes and not run.splices and not run.failures:
        return
    facts, evidence = [], []
    if sizes:
        designs = {f"{treatment} vs {control}" for treatment, control in sizes.values()}
        if len(designs) == 1:
            facts.append(say(section, "design", designs=designs.pop()))
        else:
            facts.append(
                say(
                    section,
                    "designs",
                    designs=", ".join(
                        f"{treatment} vs {control} ({name})"
                        for name, (treatment, control) in sorted(sizes.items())
                    ),
                )
            )
    if "LeafCutter" in run.ran and ARGS["leafcutter_below_defaults"]:
        facts.append(say(section, "leafcutter"))
    if run.splices:
        values = pd.Series(run.splices)
        spread = splices_text(values.min())
        if values.min() != values.max():
            spread += f" to {splices_text(values.max())}"
        facts.append(
            say(section, "splices", median=splices_text(values.median()), range=spread)
        )
    smallest = min((size for pair in sizes.values() for size in pair), default=None)
    if smallest is not None and run.ran & DTU and smallest <= 2:
        evidence.append(say(section, "evidence_dtu"))
    if smallest is not None and "LeafCutter" in run.ran and smallest < 4:
        evidence.append(say(section, "evidence_leafcutter"))
    text = " ".join(facts)
    if evidence:
        sentence = "; ".join(evidence)
        text += f" {sentence[0].upper()}{sentence[1:]}. " + say(section, "exploratory")
    content = say(section, "facts", message=text) if text else ""
    for failure in run.failures:
        content += say(
            section,
            "failure",
            **{key: html.escape(value) for key, value in failure._asdict().items()},
        )
    emit(section, content)


def design_problems(section, run, n_treatment, n_control):
    """What a contrast cannot run, where it runs below documented guidance, and the evidence that applies to it."""
    smallest, pair, sizes = (
        min(n_treatment, n_control),
        (n_treatment, n_control),
        {"n_treatment": n_treatment, "n_control": n_control},
    )
    cannot, below, evidence = [], [], []
    variance = [
        name
        for analysis, name in CONFIG["variance_analyses"].items()
        if analysis in run.ran
    ]
    if smallest < 2 and variance:
        cannot.append(say(section, "variance", names=", ".join(variance), **sizes))
    for analysis, option, key in (
        ("SUPPA2", "suppa_paired", "suppa_paired"),
        ("rMATS", "rmats_paired_stats", "rmats_paired"),
    ):
        if analysis in run.ran and ARGS[option] and n_treatment != n_control:
            cannot.append(say(section, key, **sizes))
    if "LeafCutter" in run.ran:
        if smallest < ARGS["leafcutter_min_samples_per_condition"]:
            cannot.append(say(section, "leafcutter_minimum", **sizes))
        if ARGS["leafcutter_below_defaults"]:
            below.append(say(section, "leafcutter_defaults"))
        if smallest < 4:
            below.append(say(section, "leafcutter_few", **sizes))
            evidence.append(say(section, "evidence_leafcutter_few"))
        if pair == (2, 2) and ARGS["leafcutter_min_samples_per_intron"] < 4:
            below.append(say(section, "leafcutter_i"))
            evidence.append(say(section, "evidence_leafcutter_i"))
    needed = max(
        ARGS[f"min_samps_{name}"]
        for name in ("gene_expr", "feature_expr", "feature_prop")
    )
    if "DRIMSeq" in run.ran and len(run.condition_of) < needed:
        cannot.append(
            say(section, "drimseq", needed=needed, total=len(run.condition_of))
        )
    if run.ran & DTU and smallest <= 3:
        evidence.append(
            say(section, "evidence_dtu_2" if smallest <= 2 else "evidence_dtu_3")
        )
    if pair == (2, 2):
        evidence.append(say(section, "evidence_permutation"))
    return cannot, below, evidence


def write_design_check(section, run, calls, counts):
    if not run.contrasts or not run.condition_of:
        return
    leafcutter, drimseq = "LeafCutter" in run.ran, "DRIMSeq" in run.ran
    shown = {
        "splices": bool(run.splices),
        "drimseq_gene": drimseq,
        "leafcutter_g": leafcutter,
        "leafcutter_i": leafcutter,
        "leafcutter_c": leafcutter,
    }
    data = {}
    for name, (n_treatment, n_control) in sorted(run.sizes().items()):
        contrast = run.contrasts[name]
        row = {
            "treatment": contrast.treatment,
            "control": contrast.control,
            "n_treatment": n_treatment,
            "n_control": n_control,
        }
        if leafcutter:
            row |= {
                f"leafcutter_{flag}": ARGS[option]
                for flag, option in CONFIG["leafcutter_flags"].items()
            }
        if drimseq:
            row["drimseq_gene"] = ARGS["min_samps_gene_expr"]
        splices = [
            run.splices[sample]
            for sample in run.samples(contrast.treatment)
            | run.samples(contrast.control)
            if sample in run.splices
        ]
        if splices:
            row["splices"] = np.median(splices) / 1e6
        cannot, below, evidence = design_problems(section, run, n_treatment, n_control)
        row["status"] = say(
            section, "cannot" if cannot else "below" if below else "meets"
        )
        details = [
            say(section, key, items="; ".join(items))
            for key, items in (("details_cannot", cannot), ("details_below", below))
            if items
        ]
        data[name] = row | {
            "details": " ".join(details),
            "evidence": "; ".join(evidence),
        }
    emit(
        section,
        data,
        headers={
            key: header
            for key, header in SECTIONS[section]["headers"].items()
            if shown.get(key, True)
        },
    )


def write_concordance(section, run, calls, counts):
    route = first_route(run.psi)
    if route is None:
        return
    psi = run.psi[route]
    informative = psi.dropna()[psi.dropna().mean(axis=1).between(0.05, 0.95)]
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
    ctx = {
        "total": f"{len(informative):,}",
        "route": routed("from", route) if len(run.psi) > 1 else "",
        "few": "",
    }
    if len(informative) < LIMITS["psi_concordance_few_events"]:
        ctx["few"] = say(section, "few")
    emit(
        section,
        correlation.where(correlation.notna(), None).to_numpy().tolist(),
        ctx,
        xcats=names,
        ycats=names,
    )


def overview(calls):
    """Tested and significant features of each analysis, route, unit and contrast, with the median effect size of the significant ones."""
    counts = (
        calls[calls["counted"]]
        .groupby(GROUP, sort=False)["significant"]
        .agg(tested="size", significant="sum")
        .reset_index()
    )
    significant = calls[calls["significant"]]
    effects = (
        significant["effect"]
        .abs()
        .groupby([significant[column] for column in GROUP])
        .median()
    )
    units = [
        UNITS[method][unit] for method, unit in zip(counts["method"], counts["unit"])
    ]
    sources = [
        counts["method"],
        counts["route"],
        [unit.get("effect_from", name) for unit, name in zip(units, counts["unit"])],
        counts["contrast"],
    ]
    counts["effect"] = effects.reindex(pd.MultiIndex.from_arrays(sources)).to_numpy()
    counts.loc[["effect_unit" not in unit for unit in units], "effect"] = np.nan
    rank = {"method": list(METHODS).index, "route": ROUTES.index}
    return counts.sort_values(
        GROUP,
        key=lambda column: (
            column.map(rank[column.name]) if column.name in rank else column
        ),
        kind="stable",
    )


def column_headers(section, row):
    """The significant, tested and effect headers of one analysis, route and unit of the cross-method table."""
    unit = UNITS[row.method][row.unit]
    ctx = ARGS | {
        "method": row.method,
        "unit": row.unit,
        "effect_unit": unit.get("effect_unit"),
        "route_text": routed("quantification", row.route),
    }
    ctx |= {
        "title": f"{row.method} {row.unit}{routed('label', row.route)}",
        "threshold": unit["threshold"].format_map(ARGS),
        "effect": unit.get("effect", CONFIG["no_effect_cutoff"]).format_map(ARGS),
        "tested_text": unit.get("tested", CONFIG["tested"]).format_map(ctx),
    }
    templates = SECTIONS[section]["headers"]
    return {
        name: {
            field: value.format_map(ctx) if field in ("title", "description") else value
            for field, value in templates[name].items()
        }
        for name in templates
    }


def write_differential(section, run, calls, counts):
    if counts.empty:
        return
    primary_routes = counts.groupby("method")["route"].agg(
        lambda routes: first_route(set(routes))
    )
    contrasts = sorted(set(counts["contrast"]))
    headers, data = {}, defaultdict(dict)
    reversed_contrasts = set(contrasts) & run.reverses.other.keys()
    if reversed_contrasts:
        headers["reciprocal"] = SECTIONS[section]["headers"]["reciprocal"]
    for contrast in reversed_contrasts:
        data[contrast]["reciprocal"] = run.reverses.other[contrast]
    for row in counts.itertuples():
        primary = row.route == primary_routes[row.method]
        key = re.sub(
            "[^a-z0-9]+", "_", f"{row.method} {row.route} {row.unit}".lower()
        ).strip("_")
        if key not in headers:
            column = column_headers(section, row)
            headers[key] = column["significant"] | {
                "hidden": not (
                    primary and METHODS[row.method].get("primary") == row.unit
                )
            }
            headers[f"{key}_tested"] = column["tested"]
            if "effect_unit" in UNITS[row.method][row.unit] and primary:
                headers[f"{key}_effect"] = column["effect"]
        data[row.contrast] |= {key: row.significant, f"{key}_tested": row.tested}
        if primary and pd.notna(row.effect):
            data[row.contrast][f"{key}_effect"] = row.effect
    ctx = {"reciprocal": say(section, "reciprocal") if reversed_contrasts else ""}
    emit(
        section,
        {contrast: data[contrast] for contrast in contrasts},
        ctx,
        headers=headers,
    )


def write_leafcutter(section, run, calls, counts):
    if run.leafcutter:
        categories = {
            key: category | {"name": category["name"].format_map(ARGS)}
            for key, category in SECTIONS[section]["categories"].items()
        }
        emit(
            section,
            {
                run.reverses.label(contrast): outcomes
                for contrast, outcomes in sorted(run.leafcutter.items())
            },
            categories=categories,
        )


def write_suppa_events(section, run, calls, counts):
    events = calls[(calls["method"] == "SUPPA2") & (calls["unit"] == "events")]
    route = first_route(set(events["route"]))
    if route is None:
        return
    ctx = {
        "route": routed("from", route) if events["route"].nunique() > 1 else "",
        "without": "",
    }
    events = events[events["route"] == route]
    found = events.groupby("contrast")["significant"].any()
    if not found.any():
        return emit(section, say(section, "none"), ctx, plot_type="html")
    if not found.all():
        ctx["without"] = say(
            section,
            "without",
            contrasts=", ".join(map(run.reverses.label, found.index[~found])),
        )

    def category(event_type, higher):
        return say(
            section,
            "category",
            event_type=event_type,
            direction="higher" if higher else "lower",
        )

    typed, minimum = (
        events[events["event_type"].notna()],
        LIMITS["suppa_type_min_significant"],
    )
    totals, shares = {}, {}
    for contrast in found.index[found]:
        of_contrast = typed[typed["contrast"] == contrast]
        significant = of_contrast[of_contrast["significant"]]
        pairs = Counter(zip(significant["event_type"], significant["effect"] > 0))
        per_type, of_type = (
            significant["event_type"].value_counts(),
            of_contrast["event_type"].value_counts(),
        )
        totals[run.reverses.label(contrast)] = {
            category(*pair): count for pair, count in pairs.items()
        }
        share = {
            category(*pair): 100 * count / of_type[pair[0]]
            for pair, count in pairs.items()
            if per_type[pair[0]] >= minimum
        }
        if share:
            shares[run.reverses.label(contrast)] = share
    colours = CONFIG["suppa_event_types"].values()
    categories = {
        category(name, higher): {"name": category(name, higher), "color": colour}
        for name, *shades in colours
        for higher, colour in zip((True, False), shades)
    }
    data, pconfig = totals, dict(SECTIONS[section]["pconfig"])
    if shares:
        data, categories = [totals, shares], [categories, categories]
        pconfig["data_labels"] = [
            {key: value.format_map(ARGS) for key, value in label.items()}
            for label in SECTIONS[section]["data_labels"]
        ]
    emit(section, data, ctx, pconfig=pconfig, categories=categories)


def ordered(methods):
    return [method for method in METHODS if method in set(methods)]


def write_agreement(section, run, calls, counts):
    primary = calls[
        calls["unit"]
        == calls["method"].map(lambda method: METHODS[method].get("primary"))
    ]
    calling = primary[~primary["contrast"].isin(run.reverses.later)]
    calling = calling[calling.groupby("contrast")["method"].transform("nunique") >= 2]
    if calling.empty:
        return
    two_routes = calling.groupby(["contrast", "method"])["route"].nunique().gt(1).any()
    hits = calling[calling["significant"]].dropna(subset=["gene"])
    hits = hits.assign(gene=hits["gene"].astype(str).str.split(",")).explode("gene")
    hits = hits.assign(gene=hits["gene"].map(run.genes.key)).dropna(subset=["gene"])
    per_method = (
        hits.groupby(["contrast", "gene", "method"])["route"].nunique().reset_index()
    )
    per_method["both"] = per_method["method"].where(per_method["route"] > 1)
    genes = (
        per_method.groupby(["contrast", "gene"])
        .agg(methods=("method", ordered), both=("both", ordered))
        .reset_index()
    )
    genes["families"] = genes["methods"].map(
        lambda methods: sorted({METHODS[method]["family"] for method in methods})
    )
    genes = genes.assign(
        name=genes["gene"].map(run.genes.name),
        n_families=genes["families"].map(len),
        n_methods=genes["methods"].map(len),
    )
    genes = genes[genes["n_families"] >= 2].sort_values(
        ["contrast", "n_families", "n_methods", "name"],
        ascending=[True, False, False, True],
        kind="stable",
    )
    sizes, maximum = genes.groupby("contrast").size(), LIMITS["agreement_max_genes"]
    table = {}
    for gene in genes.groupby("contrast").head(maximum).itertuples():
        table[f"{gene.name} ({gene.contrast})"] = {
            "gene_id": gene.gene,
            "n_families": gene.n_families,
            "n_methods": gene.n_methods,
            "families": ", ".join(gene.families),
            "methods": ", ".join(gene.methods),
        }
        if two_routes:
            table[f"{gene.name} ({gene.contrast})"]["both_routes"] = (
                ", ".join(gene.both) or "No"
            )
    contrasts = set(primary["contrast"])
    reverse = [
        say(
            section,
            "reverse_item",
            contrast=contrast,
            other=run.reverses.other[contrast],
        )
        for contrast in sorted(run.reverses.later & contrasts)
    ]
    capped = [
        f"{contrast} ({size:,})" for contrast, size in sizes.items() if size > maximum
    ]
    ctx = {
        "reverse": say(section, "reverse", contrasts=", ".join(reverse))
        if reverse
        else "",
        "capped": "",
    }
    if capped:
        ctx["capped"] = say(section, "capped", contrasts=", ".join(capped))
    if not table:
        return emit(section, say(section, "none"), ctx, plot_type="html")
    emit(
        section,
        table,
        ctx,
        headers={
            key: header
            for key, header in SECTIONS[section]["headers"].items()
            if key != "both_routes" or two_routes
        },
    )


def write_pvalues(section, run, calls, counts):
    bins = LIMITS["p_bins"]
    keys = list(zip(calls["method"], calls["unit"], calls["route"]))
    labels = {
        key: UNITS[key[0]][key[1]]["histogram"] + routed("label", key[2])
        for key in set(keys)
        if "histogram" in UNITS[key[0]][key[1]]
    }
    pvalues = calls.assign(series=[labels.get(key) for key in keys])
    pvalues = pvalues[pvalues["series"].notna() & pvalues["p"].between(0, 1)]
    pvalues = pvalues[
        pvalues.groupby(["contrast", "series"], sort=False)["p"].transform("size")
        >= LIMITS["p_min_features"]
    ]
    contrasts = sorted(set(pvalues["contrast"]))
    if not contrasts:
        return

    def histogram(p):
        binned = np.bincount(
            np.minimum((p * bins).astype(int), bins - 1), minlength=bins
        )
        return {
            f"{(position + 0.5) / bins:.3f}": 100 * count / len(p)
            for position, count in enumerate(binned)
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
    pconfig = SECTIONS[section]["pconfig"]
    pconfig = pconfig | {
        "data_labels": [
            {"name": run.reverses.label(contrast), "ylab": pconfig["ylab"]}
            for contrast in contrasts
        ]
    }
    emit(section, data, {"per_bin": f"{100 / bins:g}"}, pconfig=pconfig)


WRITERS = {
    "splicing_design_check": write_design_check,
    "suppa_psi_concordance": write_concordance,
    "differential_splicing": write_differential,
    "leafcutter_cluster_test_outcomes": write_leafcutter,
    "suppa_significant_events": write_suppa_events,
    "splicing_method_agreement": write_agreement,
    "splicing_pvalue_histograms": write_pvalues,
}


def main():
    run = Run()
    for label, path in sorted(
        zip(LABELS, STAGED),
        key=lambda item: (
            PRELOAD.index(item[0]["kind"])
            if item[0]["kind"] in PRELOAD
            else len(PRELOAD)
        ),
    ):
        kind = KINDS.get(label["kind"], {})
        if "analysis" in kind:
            run.ran.add(kind["analysis"])
        with fail_soft(run, kind.get("name", label["kind"]), os.path.basename(path)):
            for frame in PARSERS[label["kind"]](run, label, path) or []:
                run.frames.append(
                    frame.assign(
                        route=label["route"],
                        contrast=frame.get("contrast", label["id"]),
                    )
                )
    psi = run.psi.get(first_route(run.psi))
    for sample in psi if psi is not None else []:
        run.general_stats[sample]["suppa_events_quantified"] = (
            100 * psi[sample].notna().sum() / len(psi)
        )
    run.reverses = Reverses(run.contrasts)
    columns = [*GROUP, "significant", "gene", "p", "effect", "event_type", "counted"]
    calls = pd.concat(
        [frame.reindex(columns=columns) for frame in run.frames]
        or [pd.DataFrame(columns=columns)],
        ignore_index=True,
    )
    calls = calls.astype(
        {"counted": bool, "significant": bool, "p": float, "effect": float}
    )
    counts = overview(calls)
    if run.general_stats:
        write_general_stats(run)
    for section, write in WRITERS.items():
        with fail_soft(run, SECTIONS[section]["section_name"], f"{section}_mqc.json"):
            write(section, run, calls, counts)
    write_verdict("splicing_verdict", run, calls, counts)
    versions = {"python": platform.python_version()} | {
        package.__name__: package.__version__ for package in (np, pd, yaml)
    }
    with open("versions.yml", "w") as handle:
        yaml.safe_dump({PROCESS: versions}, handle)


if __name__ == "__main__":
    main()
