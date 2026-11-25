#!/usr/bin/env python

# Summarise the per-event, per-intron, per-exon-bin and per-transcript tables of the
# splicing analyses into small MultiQC custom-content files: per-sample columns for the
# General Statistics table, and one file per report section, all grouped under a Splicing
# parent section.

import argparse
import csv
import gzip
import json
import math
import operator
import re
import shlex
import statistics
import textwrap
from collections import defaultdict

labels = json.loads(r"""${groovy.json.JsonOutput.toJson(labels)}""")
staged = json.loads(r"""${groovy.json.JsonOutput.toJson((results instanceof List ? results : [results]).collect { result -> result.toString() })}""")

PARENT = {
    "parent_id": "splicing",
    "parent_name": "Splicing",
    "parent_description": "Whether the data supported the splicing analyses, what they found and whether they agree, computed by the pipeline from their result tables.",
}

PERCENT = {"suffix": "%", "min": 0, "max": 100, "scale": "RdYlGn"}
GENERAL_STATS_HEADERS = {
    "star_annotated_splices": {
        **PERCENT,
        "title": "Annotated Splices",
        "namespace": "STAR",
        "description": "Percentage of the splices STAR found that are annotated junctions (Number of splices: Annotated (sjdb) over Total, from Log.final.out). Below about 95% suggests a genome/annotation mismatch.",
    },
    "star_splices": {
        "title": "Spliced Reads",
        "namespace": "STAR",
        "description": "Spliced reads per sample from STAR, in millions. No validated minimum exists; in one benchmark, detection of validated events was mostly stable from 40-50 M total reads per sample (Mehmood et al. 2020). See Junction Saturation for whether this sample's junction detection has levelled off.",
        "suffix": " M",
        "format": "{:,.2f}",
        "min": 0,
        "scale": "Blues",
    },
    "dexseq_reads_in_bins": {
        **PERCENT,
        "title": "Reads in Exon Bins",
        "namespace": "DEXSeq",
        "description": "Approximate percentage of the reads dexseq_count counted into exon bins, against those it counted as _empty, _ambiguous, _lowaqual or _notaligned; a read overlapping several exon bins of one gene counts once for each. Low values usually mean dexseq_count ran with the wrong strandedness or a mismatched annotation.",
    },
    "suppa_events_quantified": {
        **PERCENT,
        "title": "Events Quantified",
        "namespace": "SUPPA2",
        "description": "Percentage of the SUPPA2 local splicing events with a PSI value in the sample, from the STAR + Salmon route when both routes ran",
        "hidden": True,
    },
    "leafcutter_junction_reads": {
        "title": "Junction Reads",
        "namespace": "LeafCutter",
        "description": "Reads on the introns of the LeafCutter clusters in the sample, a read spanning several introns counting once for each. leafcutter_ds tests a cluster only in samples with at least --leafcutter_ds_min_coverage (-c) reads on it",
        "format": "{:,.0f}",
        "min": 0,
        "scale": "Blues",
    },
    "leafcutter_introns_detected": {
        **PERCENT,
        "title": "Introns Detected",
        "namespace": "LeafCutter",
        "description": "Percentage of the introns in the LeafCutter clusters with at least one read in the sample",
        "hidden": True,
    },
    "leafcutter_reads_strandless": {
        "title": "Reads Strandless",
        "namespace": "LeafCutter",
        "description": "Percentage of the junction reads of an unstranded sample on junctions that got no strand from the aligner XS tag, the splice motif or the annotation. LeafCutter discards these junctions without a message, so more than about 5% is worth checking",
        "suffix": "%",
        "min": 0,
        "max": 100,
        "scale": "OrRd",
        "hidden": True,
        "cond_formatting_rules": {"warn": [{"gt": 5}]},
    },
}

METHOD_ORDER = [
    "DEXSeq DEU",
    "edgeR DEU",
    "DEXSeq DTU",
    "DEXSeq DTU + stageR",
    "SUPPA2",
    "IsoformSwitchAnalyzeR",
    "LeafCutter",
    "ASpli",
    "rMATS",
]

# The column of each analysis shown in the cross-method table, from its primary route
PRIMARY_UNITS = {
    "DEXSeq DEU": "genes",
    "edgeR DEU": "genes",
    "DEXSeq DTU + stageR": "genes",
    "SUPPA2": "events",
    "IsoformSwitchAnalyzeR": "genes with switches",
    "LeafCutter": "clusters",
    "ASpli": "bins",
    "rMATS": "JC events",
}

# Analyses that measure the same thing from the same quantification are one line of evidence
FAMILIES = {
    "DEXSeq DEU": "Exon usage",
    "edgeR DEU": "Exon usage",
    "DEXSeq DTU + stageR": "Transcript usage",
    "IsoformSwitchAnalyzeR": "Transcript usage",
    "SUPPA2": "SUPPA2 events",
    "LeafCutter": "LeafCutter clusters",
    "ASpli": "ASpli bins",
    "rMATS": "rMATS events",
}

ROUTE_ORDER = ["STAR", "Salmon", ""]
ROUTE_LABELS = {"STAR": "STAR + Salmon", "Salmon": "Salmon", "": ""}

NO_EFFECT_CUTOFF = "no effect-size cutoff"
# The FDR cutoff of the rMATS summary.txt, which the pipeline does not change
RMATS_FDR = 0.05
P_BINS = 20
P_MIN_FEATURES = 100
LEAFCUTTER_DEFAULTS = {"g": 3, "i": 5, "c": 20}
DTU_ANALYSES = {"DRIMSeq", "IsoformSwitchAnalyzeR"}
KIND_ANALYSES = {
    "dexseq_deu_bins": "DEXSeq DEU",
    "dexseq_deu_genes": "DEXSeq DEU",
    "edger_exons": "edgeR DEU",
    "edger_genes": "edgeR DEU",
    "dexseq_dtu_transcripts": "DRIMSeq",
    "dexseq_dtu_genes": "DRIMSeq",
    "stager": "DRIMSeq",
    "suppa_events": "SUPPA2",
    "suppa_isoforms": "SUPPA2",
    "isa_features": "IsoformSwitchAnalyzeR",
    "leafcutter_numers": "LeafCutter",
    "leafcutter_ds": "LeafCutter",
    "aspli_bins": "ASpli",
    "aspli_clusters": "ASpli",
    "rmats_summary": "rMATS",
    "rmats_jc": "rMATS",
}
PSI_CONCORDANCE_MIN_EVENTS = 50
PSI_CONCORDANCE_FEW_EVENTS = 200
# Pure Python correlations are quadratic in samples, so the concordance uses a sample of events
PSI_CONCORDANCE_MAX_EVENTS = 5000
SUPPA_TYPE_MIN_SIGNIFICANT = 10
AGREEMENT_MAX_GENES = 100

DEXSEQ_COUNTERS = {
    "_empty": "empty",
    "_ambiguous": "ambiguous",
    "_ambiguous_readpair_position": "ambiguous_readpair_position",
    "_lowaqual": "lowaqual",
    "_notaligned": "notaligned",
}

SUPPA_EVENT_TYPES = {
    "SE": "Skipped exon",
    "A5": "Alternative 5' splice site",
    "A3": "Alternative 3' splice site",
    "MX": "Mutually exclusive exons",
    "RI": "Retained intron",
    "AF": "Alternative first exon",
    "AL": "Alternative last exon",
}
SUPPA_EVENT_TYPE_COLOURS = {
    "Skipped exon": "#437bb1",
    "Alternative 5' splice site": "#e8853a",
    "Alternative 3' splice site": "#3a9a5b",
    "Mutually exclusive exons": "#8e5ea2",
    "Retained intron": "#c0392b",
    "Alternative first exon": "#2a9d8f",
    "Alternative last exon": "#9c6b30",
}
# Local event IDs are `<gene>;<type>:<coordinates>`, isoform IDs `<gene>;<transcript>`
SUPPA_EVENT_ID = re.compile(r"^[^;\\t]+;(" + "|".join(SUPPA_EVENT_TYPES) + r"):\\S+")

# The status leafcutter_ds gives a cluster it did not test, with a key and label for each
LEAFCUTTER_SKIP_STATUSES = {
    "Too many introns in cluster": ("too_many_introns", "Not tested: too many introns"),
    "<=1 junction in cluster": ("one_junction", "Not tested: one junction"),
    "<=1 sample with coverage>0": ("one_sample_covered", "Not tested: one sample with reads"),
    "<=1 sample with coverage>min_coverage": ("one_sample_min_coverage", "Not tested: one sample above min coverage"),
    "<2 introns used in >=min_samples_per_intron samples": ("too_few_introns_used", "Not tested: too few introns in enough samples"),
    "Not enough valid samples": ("not_enough_samples", "Not tested: not enough valid samples"),
}
LEAFCUTTER_NOT_TESTED_COLOURS = ["#e0e0e0", "#cccccc", "#b1b1b1", "#969696", "#7a7a7a", "#5e5e5e", "#424242", "#262626"]


def parse_ext_args(args_string):
    if args_string in (None, "", "null"):
        args_string = ""

    def flag(value):
        return str(value).lower() == "true"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fdr", type=float, default=0.05, help="Adjusted p-value cutoff of the tools without one of their own")
    parser.add_argument("--stager_alpha", type=float, default=0.05, help="Target OFDR of stageR, --stager_alpha")
    parser.add_argument("--aspli_bin_fdr", type=float, default=0.05, help="Bin FDR cutoff of ASpli, --aspli_bin_fdr")
    parser.add_argument("--suppa_alpha", type=float, default=0.05, help="p-value cutoff of SUPPA2 diffSplice")
    parser.add_argument("--suppa_dpsi", type=float, default=0.05, help="Absolute dPSI cutoff of SUPPA2, the --clusterevents_dpsithreshold of clusterEvents")
    parser.add_argument("--leafcutter_min_samples_per_condition", type=int, default=5, help="Samples each condition needs for the pipeline to test a LeafCutter contrast")
    parser.add_argument("--suppa_gene_correction", type=flag, default=True, help="Whether diffSplice corrected its p-values by gene")
    parser.add_argument("--suppa_paired", type=flag, default=False, help="Whether diffSplice ran a paired test")
    parser.add_argument("--isa_alpha", type=float, default=0.05, help="Switch q-value cutoff of IsoformSwitchAnalyzeR")
    parser.add_argument("--isa_dif", type=float, default=0.1, help="Absolute dIF cutoff of IsoformSwitchAnalyzeR")
    parser.add_argument("--leafcutter_min_samples_per_group", type=int, default=3)
    parser.add_argument("--leafcutter_min_samples_per_intron", type=int, default=5)
    parser.add_argument("--leafcutter_min_coverage", type=int, default=20)
    parser.add_argument("--min_samps_gene_expr", type=int, default=4)
    parser.add_argument("--min_samps_feature_expr", type=int, default=2)
    parser.add_argument("--min_samps_feature_prop", type=int, default=2)
    parser.add_argument("--rmats_paired_stats", type=flag, default=False)
    return parser.parse_args(shlex.split(args_string))


def to_float(value):
    """A number, or None for the missing values R and Python write."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def fmt(number):
    return f"{number:g}"


def percentage(part, whole):
    return 100 * part / whole if whole else None


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, newline="")


def read_delimited(path, delimiter="\\t"):
    """Rows of a table with a header, as dictionaries, one at a time."""
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def route_label(route):
    return ROUTE_LABELS.get(route, route)


def with_route(name, route):
    return f"{name} ({route_label(route)})" if route else name


def first_route(routes):
    """The STAR and Salmon routes quantify the same samples, so per-sample summaries show one."""
    for route in ROUTE_ORDER:
        if route in routes:
            return route
    return None


def write_section(section_id, name, description, plot_type, data, helptext, pconfig=None, **extra):
    content = {
        "id": section_id,
        **PARENT,
        "section_name": name,
        "description": description,
        "helptext": textwrap.dedent(helptext).strip(),
        "plot_type": plot_type,
    }
    if pconfig:
        content["pconfig"] = pconfig
    content.update(extra)
    content["data"] = data
    with open(f"{section_id}_mqc.json", "w") as handle:
        json.dump(content, handle, indent=4)


def splices_text(count):
    return f"{count / 1e6:.2f} M" if count >= 1e5 else f"{count:,.0f}"


class Reverses:
    """Contrasts that test the same comparison as another in reverse; the later of each pair is the reverse."""

    def __init__(self, contrasts):
        by_pair = {(contrast["treatment"], contrast["control"]): name for name, contrast in contrasts.items()}
        self.other = {
            name: by_pair[(contrast["control"], contrast["treatment"])]
            for name, contrast in contrasts.items()
            if (contrast["control"], contrast["treatment"]) in by_pair
        }

    def is_reverse(self, contrast):
        return contrast in self.other and contrast > self.other[contrast]

    def label(self, contrast):
        return f"{contrast} (reverse of {self.other[contrast]})" if self.is_reverse(contrast) else contrast


#
# Genes
#


class Genes:
    """Gene IDs without version suffixes, with the names some tools report mapped to their ID."""

    def __init__(self):
        self.names = {}
        self.ids_by_name = {}

    def load(self, path):
        ids_by_name = defaultdict(set)
        for row in read_delimited(path):
            gene_id = self.strip_version(row["gene_id"])
            self.names[gene_id] = row.get("gene_name") or gene_id
            ids_by_name[self.names[gene_id]].add(gene_id)
        self.ids_by_name = {name: ids.pop() for name, ids in ids_by_name.items() if len(ids) == 1}

    def strip_version(self, gene):
        stem, dot, suffix = gene.rpartition(".")
        return stem if dot and suffix.isdigit() and (stem.startswith("ENS") or stem in self.names) else gene

    def key(self, gene):
        gene = gene.strip().strip('"')
        if not gene or gene in ("NA", "."):
            return None
        gene = self.strip_version(gene)
        return gene if gene in self.names else self.ids_by_name.get(gene, gene)

    def name(self, gene):
        return self.names.get(gene, gene)


class Calls:
    """Tested and significant features of one table, the genes behind them, their p-values and the effect sizes of the significant ones."""

    def __init__(self):
        self.tested = 0
        self.significant = 0
        self.tested_genes = set()
        self.significant_genes = set()
        self.p_bins = [0] * P_BINS
        self.p_count = 0
        self.effects = []

    def add(self, significant, gene=None, pvalue=None, effect=None):
        self.tested += 1
        if significant:
            self.significant += 1
            if effect is not None:
                self.effects.append(abs(effect))
        self.add_gene(significant, gene)
        if pvalue is not None and 0 <= pvalue <= 1:
            self.p_bins[min(int(pvalue * P_BINS), P_BINS - 1)] += 1
            self.p_count += 1

    def add_gene(self, significant, gene):
        if gene:
            self.tested_genes.add(gene)
            if significant:
                self.significant_genes.add(gene)

    def counts(self):
        return self.tested, self.significant


def scan(path, q_column, cutoff, genes, gene_column=None, p_column=None, effect_column=None, delimiter="\\t"):
    """Features with a q-value are tested, and significant below the cutoff."""
    calls = Calls()
    for row in read_delimited(path, delimiter):
        qvalue = to_float(row.get(q_column))
        if qvalue is None:
            continue
        if effect_column and effect_column not in row:
            # DEXSeq names its fold change column after the conditions, `log2fold_<treatment>_<control>`
            effect_column = next((column for column in row if column.startswith(effect_column)), None)
        calls.add(
            qvalue < cutoff,
            genes.key(row[gene_column]) if gene_column else None,
            to_float(row.get(p_column)) if p_column else None,
            to_float(row.get(effect_column)) if effect_column else None,
        )
    return calls


#
# Per-sample tables
#


def dexseq_reads_in_bins(path):
    """Percentage of the reads counted into exon bins rather than by the counters of unassigned reads."""
    in_bins = unassigned = 0
    with open(path) as handle:
        for line in handle:
            feature, _, count = line.rstrip("\\n").partition("\\t")
            if feature in DEXSEQ_COUNTERS:
                unassigned += int(count)
            elif feature and not feature.startswith("_"):
                in_bins += int(count)
    return percentage(in_bins, in_bins + unassigned)


def star_splices(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            name, _, value = line.partition("|")
            values[name.strip()] = value.strip()
    return to_float(values.get("Number of splices: Total")), to_float(values.get("Number of splices: Annotated (sjdb)"))


def leafcutter_introns(path):
    """Clustered introns with a read, and junction reads on clustered introns, of each sample."""
    with open_text(path) as handle:
        samples = handle.readline().split()
        detected = [0] * len(samples)
        reads = [0] * len(samples)
        introns = 0
        for line in handle:
            counts = line.split()[1:]
            if not counts:
                continue
            introns += 1
            for index, count in enumerate(counts):
                count = int(count)
                reads[index] += count
                if count > 0:
                    detected[index] += 1
    return {sample: {"introns_detected": percentage(detected[index], introns), "reads": reads[index]} for index, sample in enumerate(samples)}


def suppa_psi(path):
    """
    Events with a PSI value in each sample, and an evenly spaced sample of the informative
    events (quantified in every sample, mean PSI 0.05 to 0.95), at most twice the
    concordance cap, with their total.
    """
    with open(path) as handle:
        samples = handle.readline().rstrip("\n").split("\t")
        quantified = [0] * len(samples)
        informative = []
        informative_total = 0
        stride = 1
        events = 0
        for line in handle:
            values = line.rstrip("\n").split("\t")[1:]
            if not values:
                continue
            events += 1
            # SUPPA2 releases before 2.0 wrote -1 for a missing PSI
            psis = [psi if psi is not None and psi >= 0 else None for psi in map(to_float, values)]
            for index, psi in enumerate(psis):
                if psi is not None:
                    quantified[index] += 1
            # Events near 0 or 1 in every sample would make every pair of samples look alike
            if None in psis or not 0.05 <= sum(psis) / len(psis) <= 0.95:
                continue
            if informative_total % stride == 0:
                informative.append(psis)
                if len(informative) >= 2 * PSI_CONCORDANCE_MAX_EVENTS:
                    informative = informative[::2]
                    stride *= 2
            informative_total += 1
    return {"samples": samples, "quantified": quantified, "events": events, "informative": informative, "informative_total": informative_total}


#
# Result tables
#


def suppa_event_type(event_id):
    match = SUPPA_EVENT_ID.match(event_id)
    return SUPPA_EVENT_TYPES[match.group(1)] if match else None


def suppa_diffsplice(path, alpha, dpsi_cutoff, genes):
    """Events of a SUPPA2 dPSI table tested (with a dPSI and a p-value) and significant, by type and direction."""
    calls = Calls()
    tested_by_type = defaultdict(int)
    significant_by_type = defaultdict(int)
    with open(path) as handle:
        handle.readline()
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 3:
                continue
            dpsi, pvalue = to_float(fields[1]), to_float(fields[2])
            if dpsi is None or pvalue is None:
                continue
            significant = pvalue < alpha and abs(dpsi) > dpsi_cutoff
            calls.add(significant, genes.key(fields[0].split(";")[0]), pvalue, dpsi)
            event_type = suppa_event_type(fields[0])
            if event_type:
                tested_by_type[event_type] += 1
                if significant:
                    significant_by_type[(event_type, dpsi > 0)] += 1
    return {"calls": calls, "tested_by_type": tested_by_type, "significant_by_type": significant_by_type}


def stager_levels(path, alpha, genes):
    """
    Genes from the screening stage, and transcripts from the confirmation stage, which
    only tests the transcripts of the genes that pass screening.
    """
    gene_calls = Calls()
    transcript_calls = Calls()
    screened = {}
    for row in read_delimited(path):
        gene_padj = to_float(row["gene"])
        if gene_padj is None:
            continue
        screened[row["geneID"]] = gene_padj < alpha
        if gene_padj < alpha:
            tx_padj = to_float(row["transcript"])
            transcript_calls.add(tx_padj is not None and tx_padj < alpha)
    for gene_id, passed in screened.items():
        gene_calls.add(passed, genes.key(gene_id))
    return gene_calls, transcript_calls


def isoformswitchanalyzer(path, contrasts, alpha, dif, genes):
    """Isoforms and genes with a switch q-value, by the contrast of their comparison."""
    by_pair = {(contrast["control"], contrast["treatment"]): name for name, contrast in contrasts.items()}
    isoforms = defaultdict(Calls)
    gene_switches = defaultdict(lambda: defaultdict(bool))
    for row in read_delimited(path, ","):
        qvalue = to_float(row["isoform_switch_q_value"])
        if qvalue is None:
            continue
        pair = (row["condition_1"], row["condition_2"])
        contrast = by_pair.get(pair, f"{pair[1]}-{pair[0]}")
        delta = to_float(row["dIF"])
        significant = qvalue < alpha and delta is not None and abs(delta) > dif
        isoforms[contrast].add(significant, effect=delta)
        gene_switches[contrast][genes.key(row["gene_id"])] |= significant
    gene_calls = {}
    for contrast, switches in gene_switches.items():
        gene_calls[contrast] = Calls()
        for gene, switched in switches.items():
            gene_calls[contrast].add(switched, gene)
        gene_calls[contrast].effects = isoforms[contrast].effects
    return isoforms, gene_calls


def leafcutter_outcomes(path, fdr, genes):
    """Clusters of a LeafCutter cluster significance table, by test outcome."""
    outcomes = defaultdict(int)
    calls = Calls()
    for row in read_delimited(path):
        status = row["status"]
        if status == "Success":
            padj = to_float(row["p.adjust"])
            if padj is None:
                outcomes["no_pvalue"] += 1
                continue
            significant = padj < fdr
            outcomes["significant" if significant else "not_significant"] += 1
            calls.add(significant, None, to_float(row.get("p")))
            for gene in (row.get("genes") or "").split(","):
                calls.add_gene(significant, genes.key(gene))
        elif status in LEAFCUTTER_SKIP_STATUSES:
            outcomes[LEAFCUTTER_SKIP_STATUSES[status][0]] += 1
        elif status == "timeout":
            outcomes["timeout"] += 1
        else:
            outcomes["error"] += 1
    return dict(outcomes), calls


def rmats_summary(path):
    tested = significant = 0
    for row in read_delimited(path):
        tested += int(row["TotalEventsJC"])
        significant += int(row["SignificantEventsJC"])
    return tested, significant


def rmats_jc(path, calls, genes):
    for row in read_delimited(path):
        fdr = to_float(row.get("FDR"))
        if fdr is not None:
            calls.add(fdr <= RMATS_FDR, genes.key(row["GeneID"]), to_float(row.get("PValue")), to_float(row.get("IncLevelDifference")))


#
# Sections
#


def write_general_stats(general_stats):
    columns = [column for column in GENERAL_STATS_HEADERS if any(values.get(column) is not None for values in general_stats.values())]
    content = {
        "id": "splicing_sample_metrics",
        "plot_type": "generalstats",
        "headers": {column: GENERAL_STATS_HEADERS[column] for column in columns},
        "data": {sample: {column: values[column] for column in columns if values.get(column) is not None} for sample, values in sorted(general_stats.items())},
    }
    with open("splicing_sample_metrics_mqc.json", "w") as handle:
        json.dump(content, handle, indent=4)


def contrast_sizes(contrasts, conditions):
    return {name: (len(conditions.get(contrast["treatment"], ())), len(conditions.get(contrast["control"], ()))) for name, contrast in contrasts.items()}


def range_text(values, formatter=str):
    low, high = min(values), max(values)
    return formatter(low) if low == high else f"{formatter(low)} to {formatter(high)}"


REFERENCES = """
    References:

    * Love MI, Soneson C, Patro R (2018). Swimming downstream: statistical analysis of
      differential transcript usage following Salmon quantification. F1000Research 7:952.
    * Mehmood A et al. (2020). Systematic evaluation of differential splicing tools for
      RNA-seq studies. Briefings in Bioinformatics 21:2052.
    * LeafCutter usage documentation:
      <https://davidaknowles.github.io/leafcutter/articles/Usage.html>
    """


def leafcutter_settings(args):
    return args.leafcutter_min_samples_per_group, args.leafcutter_min_samples_per_intron, args.leafcutter_min_coverage


def leafcutter_below_defaults(args):
    g, i, c = leafcutter_settings(args)
    return g < LEAFCUTTER_DEFAULTS["g"] or i < LEAFCUTTER_DEFAULTS["i"] or c < LEAFCUTTER_DEFAULTS["c"]


def leafcutter_flags(args):
    g, i, c = leafcutter_settings(args)
    return f"-g {g} -i {i} -c {c}"


LEAFCUTTER_DEFAULT_FLAGS = f"-g {LEAFCUTTER_DEFAULTS['g']} -i {LEAFCUTTER_DEFAULTS['i']} -c {LEAFCUTTER_DEFAULTS['c']}"


def write_verdict_section(contrasts, conditions, splices, ran, args):
    """The facts behind the differential calls, with the published evidence that applies to them."""
    sizes = contrast_sizes(contrasts, conditions)
    if not sizes and not splices:
        return
    facts, evidence = [], []
    if sizes:
        designs = {f"{n_treatment} vs {n_control}" for n_treatment, n_control in sizes.values()}
        if len(designs) == 1:
            facts.append(f"Design: {designs.pop()} samples per contrast.")
        else:
            facts.append("Design: " + ", ".join(f"{n_treatment} vs {n_control} ({name})" for name, (n_treatment, n_control) in sorted(sizes.items())) + " samples.")
    if "LeafCutter" in ran and leafcutter_below_defaults(args):
        facts.append(f"LeafCutter ran with <code>{leafcutter_flags(args)}</code>; its defaults are <code>{LEAFCUTTER_DEFAULT_FLAGS}</code>.")
    if splices:
        values = list(splices.values())
        facts.append(f"Spliced reads per sample: median {splices_text(statistics.median(values))} (range {range_text(values, splices_text)}).")
    smallest = min((size for pair in sizes.values() for size in pair), default=None)
    if smallest is not None and ran & DTU_ANALYSES and smallest <= 2:
        evidence.append("in simulations, most transcript-usage methods did not control the false discovery rate at 2 samples per group (Love et al. 2018)")
    if smallest is not None and "LeafCutter" in ran and smallest < 4:
        evidence.append("LeafCutter is documented as well calibrated down to 4 vs 4")
    message = " ".join(facts)
    if evidence:
        sentence = "; ".join(evidence)
        message += f" {sentence[0].upper()}{sentence[1:]}. Treat the calls below as exploratory."
    write_section(
        "splicing_verdict",
        "Design and Depth",
        "The replicates, analysis settings and sequencing depth behind the differential splicing calls, with the published evidence that applies to them.",
        "html",
        f'<div class="alert alert-info">{message}</div>',
        """
            The line gives the samples per contrast, any analysis settings below the tool
            defaults, and the spliced reads of each sample (STAR `Number of splices: Total`, which
            counts a read once for each junction it spans). Where a published evaluation or the
            tool documentation covers a design like this one, the line adds what it found, and
            then suggests treating the calls as exploratory.

            No published minimum number of replicates or depth applies to all splicing analyses.
            Power and error control improve with replicates (Mehmood et al. 2020). The Design
            Check below gives the details for each contrast and analysis, and the RSeQC Junction
            Saturation section shows whether junction detection has levelled off at the depth of
            each sample.
            """
        + REFERENCES,
    )


def write_design_section(contrasts, conditions, splices, ran, args):
    """Samples per condition of each contrast, against the requirements and documented guidance of the analyses that ran."""
    if not contrasts or not conditions:
        return
    total = len({sample for samples in conditions.values() for sample in samples})
    setting = {"format": "{:,.0f}", "scale": False, "hidden": True}
    headers = {
        "treatment": {"title": "Treatment", "description": "Treatment condition of the contrast", "scale": False},
        "control": {"title": "Control", "description": "Control condition of the contrast", "scale": False},
        "n_treatment": {"title": "Treatment samples", "description": "Samples of the treatment condition", "format": "{:,.0f}", "scale": False},
        "n_control": {"title": "Control samples", "description": "Samples of the control condition", "format": "{:,.0f}", "scale": False},
    }
    if splices:
        headers["splices"] = {
            "title": "Median spliced reads (M)",
            "description": "Median STAR splices (Number of splices: Total) of the samples of the contrast, in millions",
            "format": "{:,.2f}",
            "scale": False,
        }
    if "LeafCutter" in ran:
        headers["leafcutter_g"] = {**setting, "title": "LeafCutter -g", "description": "--leafcutter_ds_min_samples_per_group"}
        headers["leafcutter_i"] = {**setting, "title": "LeafCutter -i", "description": "--leafcutter_ds_min_samples_per_intron"}
        headers["leafcutter_c"] = {**setting, "title": "LeafCutter -c", "description": "--leafcutter_ds_min_coverage"}
    if "DRIMSeq" in ran:
        headers["drimseq_gene"] = {**setting, "title": "DRIMSeq min samples", "description": "--min_samps_gene_expr, applied to all the samples of the run"}
    headers["status"] = {
        "title": "Status",
        "description": "Cannot run, Runs below documented guidance, or Meets tool requirements. Meeting the requirements of the tools does not mean the contrast is adequately powered",
        "scale": False,
        "cond_formatting_rules": {"fail": [{"s_eq": "Cannot run"}], "warn": [{"s_eq": "Runs below documented guidance"}]},
    }
    headers["details"] = {"title": "Details", "description": "What the contrast is short of, and what to change", "scale": False}
    headers["evidence"] = {"title": "Evidence", "description": "Published evidence or tool documentation that applies to this design", "scale": False}

    g, i, c = leafcutter_settings(args)
    variance_analyses = [analysis for analysis in ("DEXSeq DEU", "edgeR DEU", "DRIMSeq", "SUPPA2") if analysis in ran]
    data = {}
    for name, (n_treatment, n_control) in sorted(contrast_sizes(contrasts, conditions).items()):
        contrast = contrasts[name]
        row = {"treatment": contrast["treatment"], "control": contrast["control"], "n_treatment": n_treatment, "n_control": n_control}
        smallest = min(n_treatment, n_control)
        cannot, below, evidence = [], [], []
        if smallest < 2 and variance_analyses:
            names = ", ".join("DEXSeq DTU" if analysis == "DRIMSeq" else analysis for analysis in variance_analyses)
            cannot.append(f"{names} need replicates to estimate variance ({n_treatment} and {n_control} samples)")
        if "SUPPA2" in ran and args.suppa_paired and n_treatment != n_control:
            cannot.append(f"the SUPPA2 paired test (--diffsplice_paired) needs equal groups ({n_treatment} and {n_control})")
        if "rMATS" in ran and args.rmats_paired_stats and n_treatment != n_control:
            cannot.append(f"rMATS paired statistics (--rmats_paired_stats) need equal groups ({n_treatment} and {n_control})")
        if "LeafCutter" in ran:
            row["leafcutter_g"], row["leafcutter_i"], row["leafcutter_c"] = g, i, c
            minimum = args.leafcutter_min_samples_per_condition
            if smallest < minimum:
                cannot.append(
                    f"LeafCutter not tested: the pipeline tests a LeafCutter contrast only with at least {minimum} samples per condition, "
                    f"the larger of -g and -i ({n_treatment} and {n_control}) (lower --leafcutter_ds_min_samples_per_group and "
                    "--leafcutter_ds_min_samples_per_intron or add replicates)"
                )
            if leafcutter_below_defaults(args):
                below.append(f"LeafCutter settings {leafcutter_flags(args)} are below its defaults {LEAFCUTTER_DEFAULT_FLAGS}")
            if smallest < 4:
                below.append(f"LeafCutter with fewer than 4 samples per group ({n_treatment} and {n_control})")
                evidence.append("LeafCutter: well calibrated down to 4 vs 4 (LeafCutter documentation)")
            if (n_treatment, n_control) == (2, 2) and i < 4:
                below.append(f"LeafCutter -i {i} is below the -i 4 (with -g 2) its documentation suggests for 2 vs 2")
                evidence.append(f"LeafCutter: -i {i} is below the documented 2 vs 2 setting -i 4")
        if "DRIMSeq" in ran:
            row["drimseq_gene"] = args.min_samps_gene_expr
            needed = max(args.min_samps_gene_expr, args.min_samps_feature_expr, args.min_samps_feature_prop)
            if total < needed:
                cannot.append(
                    f"the DRIMSeq filter needs {needed} samples and the run has {total} (lower --min_samps_gene_expr, --min_samps_feature_expr and --min_samps_feature_prop, or add samples)"
                )
        if ran & DTU_ANALYSES:
            if smallest <= 2:
                evidence.append("DTU: most methods did not control gene-level FDR at 2 per group (Love et al. 2018)")
            elif smallest == 3:
                evidence.append("DTU: DEXSeq + stageR OFDR about 10% vs 5% target at 3 per group (Love et al. 2018)")
        if (n_treatment, n_control) == (2, 2):
            evidence.append("Permutation checks impossible (3 distinct relabellings, minimum p = 1/3)")
        samples = conditions.get(contrast["treatment"], set()) | conditions.get(contrast["control"], set())
        contrast_splices = [splices[sample] for sample in samples if sample in splices]
        if contrast_splices:
            row["splices"] = statistics.median(contrast_splices) / 1e6
        if cannot:
            row["status"] = "Cannot run"
        elif below:
            row["status"] = "Runs below documented guidance"
        else:
            row["status"] = "Meets tool requirements"
        details = []
        if cannot:
            details.append(f"Cannot run: {'; '.join(cannot)}.")
        if below:
            details.append(f"Runs below documented guidance: {'; '.join(below)}.")
        row["details"] = " ".join(details)
        row["evidence"] = "; ".join(evidence)
        data[name] = row

    write_section(
        "splicing_design_check",
        "Design Check",
        "Whether each contrast meets the requirements of the splicing analyses that ran and the guidance in their documentation, "
        "with the published evidence that applies to its design. Meeting the requirements of the tools does not mean the contrast is adequately powered.",
        "table",
        data,
        """
            Each row is one contrast of the contrastsheet, with the samples of its treatment and
            control conditions and their median STAR spliced reads. The status has three tiers:

            * **Cannot run**: an analysis cannot test this contrast as configured. A condition
              with one sample gives DEXSeq, edgeR, DEXSeq DTU and SUPPA2 no estimate of the
              variance between replicates; the paired tests of SUPPA2 (`--diffsplice_paired`)
              and rMATS (`--rmats_paired_stats`) need as many samples in each condition; and
              LeafCutter and the DRIMSeq filter need the samples their settings ask for.
            * **Runs below documented guidance**: the analysis runs, but with settings below
              the tool defaults or with fewer samples than its documentation covers. LeafCutter
              defaults to `-g 3 -i 5 -c 20`, is documented as well calibrated down to 4 vs 4, and
              suggests `-i 4 -g 2` for 2 vs 2.
            * **Meets tool requirements**: none of the above. This is not a statement that the
              contrast is adequately powered.

            The Details column says what to change, and the Evidence column gives the published
            evaluations or documentation that apply to the design, where they exist. With 2 vs 2
            samples there are only 3 distinct ways to relabel the samples, so permutation checks
            cannot give a p-value below 1/3.

            No published minimum number of replicates applies to all splicing analyses. Power and
            error control improve with replicates (Mehmood et al. 2020); the entries above give
            the evidence that exists for specific tools.
            """
        + REFERENCES,
        pconfig={"id": "splicing_design_check_table", "title": "Splicing: Design Check", "namespace": "Design check", "col1_header": "Contrast", "sort_rows": False},
        headers=headers,
    )


def centred_ranks(values):
    """Average ranks of the values, minus their mean."""
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start
        while end + 1 < len(order) and values[order[end + 1]] == values[order[start]]:
            end += 1
        for position in range(start, end + 1):
            ranks[order[position]] = (start + end) / 2
        start = end + 1
    mean = (len(values) - 1) / 2
    return [rank - mean for rank in ranks]


def write_concordance_section(psi_by_route, condition_of):
    route = first_route(psi_by_route)
    if route is None:
        return
    psi = psi_by_route[route]
    total = psi["informative_total"]
    if len(psi["samples"]) < 2 or total < PSI_CONCORDANCE_MIN_EVENTS:
        return
    step = max(1, math.ceil(len(psi["informative"]) / PSI_CONCORDANCE_MAX_EVENTS))
    used = psi["informative"][::step]
    order = sorted(range(len(psi["samples"])), key=lambda index: (condition_of.get(psi["samples"][index], ""), psi["samples"][index]))
    names = [f"{psi['samples'][index]} ({condition_of[psi['samples'][index]]})" if psi["samples"][index] in condition_of else psi["samples"][index] for index in order]
    ranks = [centred_ranks([values[index] for values in used]) for index in order]
    norms = [math.sqrt(sum(map(operator.mul, rank, rank))) for rank in ranks]
    matrix = [
        [round(sum(map(operator.mul, ranks[i], ranks[j])) / (norms[i] * norms[j]), 3) if norms[i] and norms[j] else None for j in range(len(ranks))]
        for i in range(len(ranks))
    ]
    sampled = f", {len(used):,} of them taken evenly across the matrix" if len(used) < total else ""
    route_text = f" From the {route_label(route)} quantification." if len(psi_by_route) > 1 else ""
    few = (
        f" With under {PSI_CONCORDANCE_FEW_EVENTS} events the correlations are noisy and no block structure is expected."
        if total < PSI_CONCORDANCE_FEW_EVENTS
        else ""
    )
    write_section(
        "suppa_psi_concordance",
        "SUPPA2: PSI Concordance",
        f"Spearman correlation of PSI between samples ({total:,} events quantified in every sample, mean PSI 0.05-0.95{sampled}).{route_text}{few}",
        "heatmap",
        matrix,
        f"""
            PSI (percent spliced in) is the share of an event's expression that comes from its
            alternative form. This heatmap correlates the PSI values of each pair of samples over
            the SUPPA2 local events quantified in every sample, leaving out events with a mean PSI
            below 0.05 or above 0.95, which are nearly constant and would make every pair look
            alike. Samples are sorted by condition.

            With many events and a clear difference between the conditions, samples of the same
            condition correlate more closely with each other than with the other condition,
            giving blocks along the diagonal. With few events (under about
            {PSI_CONCORDANCE_FEW_EVENTS}) or a subtle contrast, no block structure is expected
            and single correlations are noisy; suspect mislabelling only when a sample's QC and
            correlations stand out consistently.

            A sample that correlates poorly with every other sample, and also stands out in the
            QC sections above, may have a library or quantification problem worth checking
            before trusting the differential results.
            """,
        pconfig={"id": "suppa_psi_concordance_plot", "title": "SUPPA2: PSI Concordance", "xlab": "Sample", "ylab": "Sample", "min": 0, "max": 1, "square": True},
        xcats=names,
        ycats=names,
    )


def write_differential_table(rows, reverses):
    """One row per contrast and one column per analysis, so the analyses of a contrast read
    across. Analyses count different features, so each column states its own in its header."""
    rows.sort(key=lambda row: (METHOD_ORDER.index(row["method"]), ROUTE_ORDER.index(row["route"]), row["unit"], row["contrast"]))
    primary_routes = {method: first_route({row["route"] for row in rows if row["method"] == method}) for method in METHOD_ORDER}
    contrasts = sorted({row["contrast"] for row in rows})
    headers = {}
    data = defaultdict(dict)
    if any(contrast in reverses.other for contrast in contrasts):
        headers["reciprocal"] = {
            "title": "Reverse of",
            "description": "The contrast that tests the same comparison in reverse, which is not an independent finding",
            "scale": False,
        }
        for contrast in contrasts:
            if contrast in reverses.other:
                data[contrast]["reciprocal"] = reverses.other[contrast]
    for row in rows:
        key = re.sub("[^a-z0-9]+", "_", " ".join(part for part in (row["method"], row["route"], row["unit"]) if part).lower())
        route = f" ({route_label(row['route'])})" if row["route"] else ""
        title = f"{row['method']} {row['unit']}{route}"
        primary = row["route"] == primary_routes[row["method"]]
        if key not in headers:
            route_text = f", from the {route_label(row['route'])} quantification" if row["route"] else ""
            headers[key] = {
                "title": title,
                "description": f"{row['method']}: {row['unit']} called significant ({row['threshold']}; {row['effect']}){route_text}",
                "format": "{:,.0f}",
                "scale": "Purples",
                "hidden": not (primary and PRIMARY_UNITS.get(row["method"]) == row["unit"]),
            }
            headers[f"{key}_tested"] = {
                "title": f"{title}: tested",
                "description": f"{row['method']}: {row['tested_text'] or row['unit'] + ' tested'}{route_text}",
                "format": "{:,.0f}",
                "scale": "Blues",
                "hidden": True,
            }
            if row["effect_unit"] and primary:
                headers[f"{key}_effect"] = {
                    "title": f"{title}: median {row['effect_unit']}",
                    "description": f"{row['method']}: median {row['effect_unit']} of the {row['unit']} called significant{route_text}",
                    "format": "{:,.2f}",
                    "scale": "Greens",
                    "hidden": True,
                }
        data[row["contrast"]][key] = row["significant"]
        data[row["contrast"]][f"{key}_tested"] = row["tested"]
        if row["effect_unit"] and primary and row["effects"]:
            data[row["contrast"]][f"{key}_effect"] = statistics.median(row["effects"])
    reciprocal_text = (
        " Contrasts in the Reverse of column test the same comparison in reverse, so they are not independent findings; "
        "for some tools, such as rMATS with paired statistics, their calls can differ slightly."
        if "reciprocal" in headers
        else ""
    )
    write_section(
        "differential_splicing",
        "Differential Splicing",
        "Features each splicing analysis called significant in each contrast, one column per analysis at gene level "
        "where the analysis has one. Counts are not comparable across columns: the analyses count different units "
        "with different tests and effect-size cutoffs. STAR + Salmon: Salmon quantification of STAR alignments; "
        f"Salmon: Salmon pseudo-alignment. An empty cell is a contrast the analysis did not test; the Design Check says why.{reciprocal_text}",
        "table",
        {contrast: data[contrast] for contrast in contrasts},
        """
            Each row is one contrast, and each visible column one analysis: the genes DEXSeq DEU,
            edgeR, DEXSeq DTU with stageR and IsoformSwitchAnalyzeR call significant, the events
            of SUPPA2 and rMATS, the intron clusters of LeafCutter and the bins of ASpli. Where an
            analysis ran on both quantification routes, the column is the STAR + Salmon route.
            The column descriptions give the threshold and effect-size cutoff of each analysis.

            More columns are hidden: the feature-level counts (exon bins, exons, transcripts,
            isoforms, junction clusters), the Salmon route, the number of features each analysis
            tested, and the median effect size of the significant calls (|log2 fold change| for
            DEXSeq DEU, edgeR and ASpli, |dPSI| for SUPPA2, |IncLevelDifference| for rMATS and
            |dIF| for IsoformSwitchAnalyzeR). Show them with the Configure Columns button.

            How to read it: a contrast with significant features in several analyses has
            splicing changes that do not depend on one method, if the run has the replicates and
            depth to support them (see Power and Design Check). Zero in an analysis that tested
            many features, in a well-replicated and deep run, suggests no large changes; it does
            not show there are none. Zero in one that tested few says little. Small median effect
            sizes with many calls point to subtle changes; large ones with few calls in a small
            run are often driven by single samples.

            rMATS counts junction-count (JC) events with FDR <= 0.05, the cutoff of the rMATS
            summary.txt, and with no |dPSI| cutoff, the rMATS default. DEXSeq tests the features left after its independent filtering, those with
            an adjusted p-value.
            """,
        pconfig={"id": "differential_splicing_table", "title": "Differential splicing", "namespace": "Differential splicing", "col1_header": "Contrast", "sort_rows": False},
        headers=headers,
    )


def write_leafcutter_outcomes_section(outcomes, args, reverses):
    categories = {
        "significant": {"name": f"Significant (adjusted p < {fmt(args.fdr)})", "color": "#7b4ea3"},
        "not_significant": {"name": "Not significant", "color": "#437bb1"},
        "no_pvalue": {"name": "Tested, no adjusted p-value", "color": "#a6c4e2"},
    }
    not_tested = [*LEAFCUTTER_SKIP_STATUSES.values(), ("timeout", "Not tested: timeout"), ("error", "Not tested: error")]
    for (key, label), colour in zip(not_tested, LEAFCUTTER_NOT_TESTED_COLOURS):
        categories[key] = {"name": label, "color": colour}
    write_section(
        "leafcutter_cluster_test_outcomes",
        "LeafCutter: Cluster Test Outcomes",
        f"Outcome of the differential splicing test of each intron cluster, one bar per contrast, with leafcutter_ds "
        f"run with -g {args.leafcutter_min_samples_per_group} (`--leafcutter_ds_min_samples_per_group`), "
        f"-i {args.leafcutter_min_samples_per_intron} (`--leafcutter_ds_min_samples_per_intron`) and "
        f"-c {args.leafcutter_min_coverage} (`--leafcutter_ds_min_coverage`). Significant: adjusted p < {fmt(args.fdr)}. "
        "Many 'Not tested' clusters usually mean too few samples per group or too little depth, not an absence of splicing changes. "
        "Very few clusters in total means too few junction reads for clustering (see LeafCutter Junction Reads in General Statistics).",
        "bargraph",
        {reverses.label(contrast): values for contrast, values in sorted(outcomes.items())},
        """
            LeafCutter groups introns that share a splice site into clusters, and `leafcutter_ds`
            tests each cluster for a change in intron usage between the conditions. It only tests
            clusters that pass its coverage filters. The reasons for not testing a cluster, as
            LeafCutter reports them in the `status` column, are:

            * **Too many introns**: more introns than `--max_cluster_size`
            * **One junction**: one junction in the cluster
            * **One sample with reads**: at most one sample with reads on the cluster
            * **One sample above min coverage**: at most one sample with `-c` reads
            * **Too few introns in enough samples**: fewer than two introns with reads in `-i` samples
            * **Not enough valid samples**: fewer than `-g` samples in a group reach `-c` reads
            * **Timeout** and **Error**: the model fit did not finish

            Good: many clusters, most of them tested, a share of them significant. Bad: most
            clusters not tested, usually for one of the coverage reasons, or very few clusters
            at all, which means the samples have too few junction reads for the clustering.
            Check the Design Check and the LeafCutter junction reads in General Statistics; add
            replicates or depth, or lower the settings above. The p-values are adjusted with the
            Benjamini-Hochberg method.
            """,
        pconfig={"id": "leafcutter_cluster_test_outcomes_plot", "title": "LeafCutter: Cluster Test Outcomes", "ylab": "Clusters", "cpswitch_counts_label": "Number of clusters"},
        categories=categories,
    )


def lighter(colour):
    channels = [int(colour[index : index + 2], 16) for index in (1, 3, 5)]
    return "#" + "".join(f"{(channel + 255) // 2:02x}" for channel in channels)


def write_suppa_events_section(diffsplice, significance, reverses):
    """diffsplice: (contrast, route) -> dPSI table summary of local events."""
    route = first_route({route for _, route in diffsplice})
    if route is None:
        return
    local = {contrast: stats for (contrast, stats_route), stats in sorted(diffsplice.items()) if stats_route == route}
    route_text = f" From the {route_label(route)} quantification." if len({route for _, route in diffsplice}) > 1 else ""
    description = f"Significant SUPPA2 events by type and direction ({significance}). With few significant events, differences between types are not interpretable.{route_text}"
    helptext = f"""
        Each bar is one contrast, split by the SUPPA2 event type and by whether the PSI is
        higher or lower in the treatment than in the control. SUPPA2 reports dPSI as the PSI of
        condition 2 minus that of condition 1, and the pipeline runs it with the control as
        condition 1, so a positive dPSI is a higher PSI in the treatment.

        * **SE**: skipped exon
        * **A5** / **A3**: alternative 5' / 3' splice site
        * **MX**: mutually exclusive exons
        * **RI**: retained intron
        * **AF** / **AL**: alternative first / last exon

        An event is significant with a `diffSplice` p-value below `--diffsplice_alpha` and an
        absolute dPSI above `--clusterevents_dpsithreshold`, the dPSI bound of SUPPA2
        `clusterEvents`, used here as the effect-size cutoff so that these counts match the
        events the pipeline clusters.

        With few significant events, differences between types or directions are not
        interpretable. The annotation defines many more events of some types than others, so
        where a type has at least {SUPPA_TYPE_MIN_SIGNIFICANT} significant events a second tab
        shows them as a percentage of the tested events of that type, which compares types more
        fairly than the counts. A strong bias towards one type or direction in a well-powered
        run can be biology (for example more intron retention in one condition), or a technical
        difference between the conditions, such as RNA degradation or library preparation,
        worth checking in the QC.
        """
    with_significant = {contrast: stats for contrast, stats in local.items() if stats["calls"].significant}
    if not with_significant:
        write_section("suppa_significant_events", "SUPPA2: Significant Events", description, "html", "<p>No contrast has any significant events.</p>", helptext)
        return
    without_significant = [reverses.label(contrast) for contrast in local if contrast not in with_significant]
    if without_significant:
        description += f" Contrasts without significant events, not shown: {', '.join(without_significant)}."
    categories = {}
    for event_type, colour in SUPPA_EVENT_TYPE_COLOURS.items():
        categories[f"{event_type}: higher in treatment"] = {"name": f"{event_type}: higher in treatment", "color": colour}
        categories[f"{event_type}: lower in treatment"] = {"name": f"{event_type}: lower in treatment", "color": lighter(colour)}

    def by_category(stats, as_percentage):
        totals = defaultdict(int)
        for (event_type, _), count in stats["significant_by_type"].items():
            totals[event_type] += count
        values = {}
        for (event_type, higher), count in stats["significant_by_type"].items():
            if as_percentage and totals[event_type] < SUPPA_TYPE_MIN_SIGNIFICANT:
                continue
            name = f"{event_type}: {'higher' if higher else 'lower'} in treatment"
            values[name] = 100 * count / stats["tested_by_type"][event_type] if as_percentage else count
        return values

    datasets = [{reverses.label(contrast): by_category(stats, False) for contrast, stats in with_significant.items()}]
    data_labels = [{"name": "Events", "ylab": "Significant events"}]
    percentages = {reverses.label(contrast): by_category(stats, True) for contrast, stats in with_significant.items()}
    percentages = {contrast: values for contrast, values in percentages.items() if values}
    if percentages:
        datasets.append(percentages)
        data_labels.append({"name": f"% of tested (types with {SUPPA_TYPE_MIN_SIGNIFICANT}+ significant)", "ylab": "% of the tested events of the type"})
    pconfig = {"id": "suppa_significant_events_plot", "title": "SUPPA2: Significant Events", "ylab": "Significant events", "stacking": "group", "cpswitch": False}
    if len(datasets) > 1:
        pconfig["data_labels"] = data_labels
    write_section(
        "suppa_significant_events",
        "SUPPA2: Significant Events",
        description,
        "bargraph",
        datasets if len(datasets) > 1 else datasets[0],
        helptext,
        pconfig=pconfig,
        categories=[categories] * len(datasets) if len(datasets) > 1 else categories,
    )


def write_agreement_section(gene_calls, genes, reverses):
    """gene_calls: contrast -> (method, route) -> Calls of genes."""
    usable = {
        contrast: methods
        for contrast, methods in gene_calls.items()
        if len({method for method, _ in methods}) >= 2 and not reverses.is_reverse(contrast)
    }
    if not usable:
        return
    two_routes = any((method, "STAR") in methods and (method, "Salmon") in methods for methods in usable.values() for method, _ in methods)
    rows = {}
    capped = []
    for contrast in sorted(usable):
        methods = usable[contrast]
        called_by = defaultdict(set)
        for (method, _), calls in methods.items():
            for gene in calls.significant_genes:
                called_by[gene].add(method)
        shared = [(gene, called, {FAMILIES.get(method, method) for method in called}) for gene, called in called_by.items()]
        shared = sorted(((gene, called, families) for gene, called, families in shared if len(families) >= 2), key=lambda item: (-len(item[2]), -len(item[1]), genes.name(item[0])))
        if len(shared) > AGREEMENT_MAX_GENES:
            capped.append(f"{contrast} ({len(shared):,})")
        for gene, called, families in shared[:AGREEMENT_MAX_GENES]:
            row = {
                "gene_id": gene,
                "n_families": len(families),
                "n_methods": len(called),
                "families": ", ".join(sorted(families)),
                "methods": ", ".join(method for method in METHOD_ORDER if method in called),
            }
            if two_routes:
                both = [
                    method
                    for method in METHOD_ORDER
                    if gene in methods.get((method, "STAR"), Calls()).significant_genes and gene in methods.get((method, "Salmon"), Calls()).significant_genes
                ]
                row["both_routes"] = ", ".join(both) or "No"
            rows[f"{genes.name(gene)} ({contrast})"] = row
    reverse_text = ", ".join(f"{contrast} is the reverse of {reverses.other[contrast]}" for contrast in sorted(reverses.other) if reverses.is_reverse(contrast) and contrast in gene_calls)
    description = (
        "Genes called significant by at least two independent lines of evidence in each contrast, sorted by the number of them. "
        "Analyses that measure the same thing from the same quantification count once: exon usage (DEXSeq DEU, edgeR), "
        "transcript usage (DEXSeq DTU + stageR, IsoformSwitchAnalyzeR), and SUPPA2, rMATS, LeafCutter and ASpli each on their own. "
        "An analysis that tests few genes calls few, so absence from this table is not disagreement."
    )
    if reverse_text:
        description += f" Reverse contrasts are not shown ({reverse_text})."
    if capped:
        description += f" Only the {AGREEMENT_MAX_GENES} genes with the most support are shown for: {', '.join(capped)}."
    helptext = """
        The analyses test different features (exon bins, events, transcripts, isoforms, intron
        clusters), so their calls are compared by gene: an analysis calls a gene significant
        when it calls any of the gene's features significant, and the two quantification
        routes of one analysis count once. Gene IDs lose their version suffix, and the gene
        names LeafCutter and IsoformSwitchAnalyzeR report are mapped to gene IDs through the
        annotation.

        Agreement between analyses that share the same quantification (DEXSeq DTU and
        IsoformSwitchAnalyzeR) is not independent support, so they count as one line of
        evidence in the Lines of evidence column. A shared gene is a candidate only if the run
        has adequate replicates and depth (see Power and Design Check); highly expressed genes
        are called first in small datasets. Where an analysis ran on both quantification
        routes, the last column lists the analyses that call the gene from both.

        Few or no shared genes is common: the analyses differ in power and in what they can
        detect. Check the tested counts in the Differential Splicing table before reading it as
        disagreement.
        """
    if not rows:
        write_section(
            "splicing_method_agreement",
            "Method Agreement: Shared Genes",
            description,
            "html",
            "<p>No gene is called significant by two or more independent lines of evidence in any contrast.</p>",
            helptext,
        )
        return
    headers = {
        "gene_id": {"title": "Gene ID", "scale": False},
        "n_families": {"title": "Lines of evidence", "description": "Independent lines of evidence that call the gene significant", "format": "{:,.0f}", "scale": "Purples"},
        "n_methods": {"title": "Analyses", "description": "Analyses that call the gene significant", "format": "{:,.0f}", "scale": "Purples"},
        "families": {"title": "Evidence", "description": "Independent lines of evidence that call the gene significant", "scale": False},
        "methods": {"title": "Called by", "description": "Analyses that call the gene significant", "scale": False},
    }
    if two_routes:
        headers["both_routes"] = {
            "title": "Called from both quantifications (STAR + Salmon and Salmon)",
            "description": "Analyses that call the gene significant from both quantification routes",
            "scale": False,
        }
    write_section(
        "splicing_method_agreement",
        "Method Agreement: Shared Genes",
        description,
        "table",
        rows,
        helptext,
        pconfig={"id": "splicing_method_agreement_table", "title": "Splicing: Shared Genes", "namespace": "Method agreement", "col1_header": "Gene (contrast)", "sort_rows": False},
        headers=headers,
    )


def write_pvalue_section(histograms, reverses):
    """histograms: contrast -> series name -> Calls."""
    datasets = {contrast: {name: calls for name, calls in series.items() if calls.p_count >= P_MIN_FEATURES} for contrast, series in histograms.items()}
    contrasts = sorted(contrast for contrast, series in datasets.items() if series)
    if not contrasts:
        return
    write_section(
        "splicing_pvalue_histograms",
        "P-value Distributions (Diagnostic)",
        f"For analyses with at least {P_MIN_FEATURES} tested features. A flat curve with a rise near 0 is expected; most "
        "features near 1 means the test had little power. One tab per contrast.",
        "linegraph",
        [
            {
                name: {f"{(position + 0.5) / P_BINS:.3f}": 100 * count / calls.p_count for position, count in enumerate(calls.p_bins)}
                for name, calls in datasets[contrast].items()
            }
            for contrast in contrasts
        ],
        f"""
            The raw p-values of the features each analysis tested, in {P_BINS} bins of 0.05, as a
            percentage of its tested features. Features without a real change have uniform
            p-values, about {100 / P_BINS:g}% per bin, and real changes add a rise near 0.

            * Flat with a rise near 0: the test behaves as expected and found changes.
            * Flat without a rise: few or no changes.
            * Most features near 1: the test had little power, often from too few replicates
              or low counts. Check the Design Check.
            * Rising towards 0 across the whole range: the p-values are too small, which can come
              from an unmodelled batch effect. Consider adding it to the design.

            SUPPA2 reports gene-corrected p-values, which are not uniform under the null, and
            IsoformSwitchAnalyzeR and stageR report only adjusted values, so they are not shown.
            """,
        pconfig={
            "id": "splicing_pvalue_histograms_plot",
            "title": "Splicing: P-value Distributions",
            "xlab": "p-value",
            "ylab": "% of tested features",
            "xmin": 0,
            "xmax": 1,
            "ymin": 0,
            "ymax": 30,
            "tt_suffix": "%",
            "data_labels": [{"name": reverses.label(contrast), "ylab": "% of tested features"} for contrast in contrasts],
        },
    )


def main():
    args = parse_ext_args("${task.ext.args ?: ''}")
    fdr = fmt(args.fdr)
    correction = "gene-corrected" if args.suppa_gene_correction else "uncorrected"
    suppa_significance = f"{correction} p < {fmt(args.suppa_alpha)}, |dPSI| > {fmt(args.suppa_dpsi)}"

    by_kind = defaultdict(list)
    for label, path in zip(labels, staged):
        by_kind[label["kind"]].append((label, path))

    genes = Genes()
    for _, path in by_kind.pop("gene_annotation", []):
        genes.load(path)
    contrasts = {}
    for _, path in by_kind.pop("contrasts", []):
        for row in read_delimited(path, ","):
            contrasts[row["contrast"]] = {"treatment": row["treatment"], "control": row["control"]}
    reverses = Reverses(contrasts)
    conditions = defaultdict(set)
    condition_of = {}
    for _, path in by_kind.pop("samplesheet", []):
        for row in read_delimited(path, ","):
            conditions[row["condition"]].add(row["sample"])
            condition_of[row["sample"]] = row["condition"]

    general_stats = defaultdict(dict)
    splices = {}
    suppa_psi_by_route = {}
    suppa_events = {}
    leafcutter_outcomes_by_contrast = {}
    gene_calls = defaultdict(dict)
    histograms = defaultdict(dict)
    rmats_calls = defaultdict(Calls)
    rows = []
    ran = set()

    def add_row(method, route, contrast, unit, counts, threshold, effect=NO_EFFECT_CUTOFF, tested_text=None, effect_unit=None, effects=()):
        tested, significant = counts
        rows.append(
            {
                "method": method,
                "route": route,
                "contrast": contrast,
                "unit": unit,
                "tested": tested,
                "significant": significant,
                "threshold": threshold,
                "effect": effect,
                "tested_text": tested_text,
                "effect_unit": effect_unit,
                "effects": list(effects),
            }
        )

    dexseq_tested = "after independent filtering, the {} with an adjusted p-value"
    log2fc = "|log2 fold change|"

    for kind, entries in by_kind.items():
        for label, path in entries:
            # The id is the sample of a per-sample table and the contrast of a result table
            route, contrast = label["route"], label["id"]
            if kind in KIND_ANALYSES:
                ran.add(KIND_ANALYSES[kind])
            if kind == "star_log":
                total, annotated = star_splices(path)
                if total:
                    splices[contrast] = total
                    general_stats[contrast]["star_annotated_splices"] = percentage(annotated or 0, total)
                    general_stats[contrast]["star_splices"] = total / 1e6
            elif kind == "dexseq_counts":
                general_stats[contrast]["dexseq_reads_in_bins"] = dexseq_reads_in_bins(path)
            elif kind == "suppa_psi":
                suppa_psi_by_route[route] = suppa_psi(path)
            elif kind == "leafcutter_numers":
                for sample, values in leafcutter_introns(path).items():
                    general_stats[sample]["leafcutter_junction_reads"] = values["reads"]
                    general_stats[sample]["leafcutter_introns_detected"] = values["introns_detected"]
            elif kind == "leafcutter_strands":
                for row in read_delimited(path):
                    general_stats[row["Sample"]]["leafcutter_reads_strandless"] = to_float(row["Reads strandless (%)"])
            elif kind == "dexseq_deu_bins":
                calls = scan(path, "padj", args.fdr, genes, gene_column="groupID", p_column="pvalue", effect_column="log2fold_")
                histograms[contrast]["DEXSeq DEU exon bins"] = calls
                add_row("DEXSeq DEU", route, contrast, "exon bins", calls.counts(), f"padj < {fdr}", tested_text=dexseq_tested.format("exon bins"), effect_unit=log2fc, effects=calls.effects)
            elif kind == "dexseq_deu_genes":
                calls = scan(path, "padj", args.fdr, genes, gene_column="groupID")
                gene_calls[contrast][("DEXSeq DEU", route)] = calls
                add_row("DEXSeq DEU", route, contrast, "genes", calls.counts(), f"gene q-value < {fdr}", tested_text=dexseq_tested.format("genes"))
            elif kind == "edger_exons":
                calls = scan(path, "FDR", args.fdr, genes, p_column="P.Value", effect_column="logFC", delimiter=",")
                histograms[contrast]["edgeR DEU exons"] = calls
                add_row("edgeR DEU", route, contrast, "exons", calls.counts(), f"FDR < {fdr}", effect_unit=log2fc, effects=calls.effects)
            elif kind == "edger_genes":
                calls = scan(path, "FDR", args.fdr, genes, gene_column="Geneid", delimiter=",")
                gene_calls[contrast][("edgeR DEU", route)] = calls
                add_row("edgeR DEU", route, contrast, "genes", calls.counts(), f"Simes FDR < {fdr}")
            elif kind == "dexseq_dtu_transcripts":
                calls = scan(path, "padj", args.fdr, genes, p_column="pvalue")
                histograms[contrast][with_route("DEXSeq DTU transcripts", route)] = calls
                add_row("DEXSeq DTU", route, contrast, "transcripts", calls.counts(), f"padj < {fdr}", tested_text=dexseq_tested.format("transcripts"))
            elif kind == "dexseq_dtu_genes":
                calls = scan(path, "padj", args.fdr, genes, gene_column="groupID")
                add_row("DEXSeq DTU", route, contrast, "genes", calls.counts(), f"gene q-value < {fdr}", tested_text=dexseq_tested.format("genes"))
            elif kind == "stager":
                gene_level, transcripts = stager_levels(path, args.stager_alpha, genes)
                gene_calls[contrast][("DEXSeq DTU + stageR", route)] = gene_level
                alpha = fmt(args.stager_alpha)
                add_row("DEXSeq DTU + stageR", route, contrast, "genes", gene_level.counts(), f"screening padj < {alpha}")
                add_row(
                    "DEXSeq DTU + stageR",
                    route,
                    contrast,
                    "transcripts",
                    transcripts.counts(),
                    f"confirmation padj < {alpha}, at a target OFDR of {alpha}",
                    tested_text="transcripts of the genes that passed screening, which the confirmation stage tests",
                )
            elif kind in ("suppa_events", "suppa_isoforms"):
                stats = suppa_diffsplice(path, args.suppa_alpha, args.suppa_dpsi, genes)
                if kind == "suppa_events":
                    suppa_events[(contrast, route)] = stats
                    gene_calls[contrast][("SUPPA2", route)] = stats["calls"]
                unit = "events" if kind == "suppa_events" else "transcripts"
                add_row(
                    "SUPPA2",
                    route,
                    contrast,
                    unit,
                    stats["calls"].counts(),
                    suppa_significance.split(",")[0],
                    effect=f"|dPSI| > {fmt(args.suppa_dpsi)}",
                    effect_unit="|dPSI|" if kind == "suppa_events" else None,
                    effects=stats["calls"].effects,
                )
            elif kind == "isa_features":
                isoforms, isa_genes = isoformswitchanalyzer(path, contrasts, args.isa_alpha, args.isa_dif, genes)
                threshold, effect = f"switch q < {fmt(args.isa_alpha)}", f"|dIF| > {fmt(args.isa_dif)}"
                for isa_contrast, calls in isoforms.items():
                    add_row("IsoformSwitchAnalyzeR", route, isa_contrast, "isoforms", calls.counts(), threshold, effect=effect)
                for isa_contrast, calls in isa_genes.items():
                    gene_calls[isa_contrast][("IsoformSwitchAnalyzeR", route)] = calls
                    add_row(
                        "IsoformSwitchAnalyzeR",
                        route,
                        isa_contrast,
                        "genes with switches",
                        calls.counts(),
                        threshold,
                        effect=effect,
                        tested_text="genes with at least one isoform with a switch q-value",
                        effect_unit="|dIF|",
                        effects=calls.effects,
                    )
            elif kind == "leafcutter_ds":
                outcomes, calls = leafcutter_outcomes(path, args.fdr, genes)
                leafcutter_outcomes_by_contrast[contrast] = outcomes
                tested = sum(outcomes.get(key, 0) for key in ("significant", "not_significant", "no_pvalue"))
                histograms[contrast]["LeafCutter clusters"] = calls
                gene_calls[contrast][("LeafCutter", route)] = calls
                add_row("LeafCutter", route, contrast, "clusters", (tested, outcomes.get("significant", 0)), f"p.adjust < {fdr}")
            elif kind == "aspli_bins":
                calls = scan(path, "bin.fdr", args.aspli_bin_fdr, genes, gene_column="locus", p_column="pvalue", effect_column="logFC")
                histograms[contrast]["ASpli bins"] = calls
                gene_calls[contrast][("ASpli", route)] = calls
                add_row("ASpli", route, contrast, "bins", calls.counts(), f"bin.fdr < {fmt(args.aspli_bin_fdr)}", effect_unit=log2fc, effects=calls.effects)
            elif kind == "aspli_clusters":
                calls = scan(path, "FDR", args.fdr, genes, p_column="pvalue")
                histograms[contrast]["ASpli junction clusters"] = calls
                add_row("ASpli", route, contrast, "junction clusters", calls.counts(), f"FDR < {fdr}")
            elif kind == "rmats_summary":
                add_row(
                    "rMATS",
                    route,
                    contrast,
                    "JC events",
                    rmats_summary(path),
                    f"FDR <= {fmt(RMATS_FDR)}",
                    effect="no |dPSI| cutoff, the rMATS default",
                    tested_text="events with junction-count (JC) statistics, the primary rMATS result",
                    effect_unit="|IncLevelDifference|",
                )
            elif kind == "rmats_jc":
                rmats_jc(path, rmats_calls[contrast], genes)
            else:
                raise ValueError(f"Unknown splicing result kind '{kind}' for {path}")

    for contrast, calls in rmats_calls.items():
        gene_calls[contrast][("rMATS", "")] = calls
        histograms[contrast]["rMATS JC events"] = calls
        for row in rows:
            if row["method"] == "rMATS" and row["contrast"] == contrast:
                row["effects"] = calls.effects

    # The DEXSeq DEU effect sizes come from its exon bins, its gene table having none
    for row in rows:
        if row["method"] == "DEXSeq DEU" and row["unit"] == "genes":
            bins = next((other for other in rows if other["method"] == "DEXSeq DEU" and other["unit"] == "exon bins" and other["contrast"] == row["contrast"]), None)
            if bins:
                row["effect_unit"], row["effects"] = f"{log2fc} of significant exon bins", bins["effects"]

    psi_route = first_route(suppa_psi_by_route)
    if psi_route is not None:
        psi = suppa_psi_by_route[psi_route]
        for index, sample in enumerate(psi["samples"]):
            general_stats[sample]["suppa_events_quantified"] = percentage(psi["quantified"][index], psi["events"])

    if general_stats:
        write_general_stats(general_stats)
    write_verdict_section(contrasts, conditions, splices, ran, args)
    write_design_section(contrasts, conditions, splices, ran, args)
    write_concordance_section(suppa_psi_by_route, condition_of)
    if rows:
        write_differential_table(rows, reverses)
    if leafcutter_outcomes_by_contrast:
        write_leafcutter_outcomes_section(leafcutter_outcomes_by_contrast, args, reverses)
    write_suppa_events_section(suppa_events, suppa_significance, reverses)
    write_agreement_section(gene_calls, genes, reverses)
    write_pvalue_section(histograms, reverses)


if __name__ == "__main__":
    main()
