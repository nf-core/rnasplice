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
import platform
import re
import shlex
import textwrap
from collections import defaultdict

labels = json.loads(r"""${groovy.json.JsonOutput.toJson(labels)}""")
staged = json.loads(r"""${groovy.json.JsonOutput.toJson((results instanceof List ? results : [results]).collect { result -> result.toString() })}""")
GENERAL_STATS_FILE = "splicing_sample_metrics_mqc.tsv"

PARENT = {
    "parent_id": "splicing",
    "parent_name": "Splicing",
    "parent_description": "Summaries of the splicing analyses, computed by the pipeline from their result tables.",
}

# Column: number format
GENERAL_STATS_COLUMNS = {
    "star_annotated_splices": "{:.2f}",
    "star_splices": "{:.4f}",
    "dexseq_reads_in_bins": "{:.2f}",
    "leafcutter_junction_reads": "{:.0f}",
    "leafcutter_introns_detected": "{:.2f}",
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

ROUTE_ORDER = ["STAR", "Salmon", ""]
ROUTE_LABELS = {"STAR": "STAR + Salmon", "Salmon": "Salmon", "": ""}

NO_EFFECT_CUTOFF = "no effect-size cutoff"
RMATS_FDR = 0.05
P_BINS = 20

DEXSEQ_COUNTERS = {
    "_empty": ("empty", "No exon bin (_empty)", "#b1b1b1"),
    "_ambiguous": ("ambiguous", "Ambiguous (_ambiguous)", "#e8853a"),
    "_ambiguous_readpair_position": ("ambiguous_readpair_position", "Ambiguous pair position", "#f2c18d"),
    "_lowaqual": ("lowaqual", "Low alignment quality (_lowaqual)", "#c0392b"),
    "_notaligned": ("notaligned", "Not aligned (_notaligned)", "#7a7a7a"),
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
# Pure Python correlations are quadratic in samples, so the concordance uses a sample of events
PSI_CONCORDANCE_MAX_EVENTS = 5000

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
AGREEMENT_MAX_GENES = 100


def parse_ext_args(args_string):
    if args_string in (None, "", "null"):
        args_string = ""

    def flag(value):
        return str(value).lower() == "true"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fdr", type=float, default=0.05, help="Adjusted p-value cutoff of the tools without one of their own")
    parser.add_argument("--stager_alpha", type=float, default=0.05, help="Target OFDR of stageR, the alpha of run_stager.R")
    parser.add_argument("--suppa_alpha", type=float, default=0.05, help="p-value cutoff of SUPPA2 diffSplice")
    parser.add_argument("--suppa_dpsi", type=float, default=0.05, help="Absolute dPSI cutoff of SUPPA2")
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


def write_section(section_id, name, description, plot_type, data, helptext=None, pconfig=None, **extra):
    content = {"id": section_id, **PARENT, "section_name": name, "description": description, "plot_type": plot_type}
    if helptext:
        content["helptext"] = textwrap.dedent(helptext).strip()
    if pconfig:
        content["pconfig"] = pconfig
    content.update(extra)
    content["data"] = data
    with open(f"{section_id}_mqc.json", "w") as handle:
        json.dump(content, handle, indent=4)


def first_route(by_route):
    """The STAR and Salmon routes quantify the same samples, so per-sample plots show one."""
    for route in ROUTE_ORDER:
        if route in by_route:
            return route
    return None


#
# Genes
#


class Genes:
    """Gene IDs without version suffixes, with the symbols some tools report mapped to their ID."""

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
        return self.names.get(gene, "")


class Calls:
    """Tested and significant features of one table, the genes behind them, and their p-values."""

    def __init__(self):
        self.tested = 0
        self.significant = 0
        self.tested_genes = set()
        self.significant_genes = set()
        self.p_bins = [0] * P_BINS
        self.p_count = 0

    def add(self, significant, gene=None, pvalue=None):
        self.tested += 1
        if significant:
            self.significant += 1
        if gene:
            self.tested_genes.add(gene)
            if significant:
                self.significant_genes.add(gene)
        if pvalue is not None and 0 <= pvalue <= 1:
            self.p_bins[min(int(pvalue * P_BINS), P_BINS - 1)] += 1
            self.p_count += 1

    def counts(self):
        return self.tested, self.significant


def scan(path, q_column, cutoff, genes, gene_column=None, p_column=None, delimiter="\\t"):
    """Features with a q-value are tested, and significant below the cutoff."""
    calls = Calls()
    for row in read_delimited(path, delimiter):
        qvalue = to_float(row.get(q_column))
        if qvalue is None:
            continue
        gene = genes.key(row[gene_column]) if gene_column else None
        calls.add(qvalue < cutoff, gene, to_float(row.get(p_column)) if p_column else None)
    return calls


#
# Per-sample tables
#


def dexseq_read_assignment(path):
    """Reads counted into exon bins and into each counter of the reads dexseq_count.py could not assign."""
    counts = {"in_bins": 0}
    with open(path) as handle:
        for line in handle:
            feature, _, count = line.rstrip("\\n").partition("\\t")
            if not feature:
                continue
            if feature in DEXSEQ_COUNTERS:
                counts[DEXSEQ_COUNTERS[feature][0]] = int(count)
            elif not feature.startswith("_"):
                counts["in_bins"] += int(count)
    return counts


def star_splices(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            name, _, value = line.partition("|")
            values[name.strip()] = value.strip()
    total = to_float(values.get("Number of splices: Total"))
    annotated = to_float(values.get("Number of splices: Annotated (sjdb)"))
    return total, annotated


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
    return {sample: {"detected": detected[index], "total": introns, "reads": reads[index]} for index, sample in enumerate(samples)}


def suppa_psi(path):
    """Events with a PSI value in each sample, and the PSI vectors of the events quantified in all of them."""
    with open(path) as handle:
        samples = handle.readline().rstrip("\\n").split("\\t")
        quantified = [0] * len(samples)
        complete = []
        events = 0
        for line in handle:
            values = line.rstrip("\\n").split("\\t")[1:]
            if not values:
                continue
            events += 1
            psis = [to_float(value) for value in values]
            for index, psi in enumerate(psis):
                # SUPPA2 releases before 2.0 wrote -1 for a missing PSI
                if psi is not None and psi >= 0:
                    quantified[index] += 1
            if all(psi is not None and psi >= 0 for psi in psis):
                complete.append(psis)
    return {"samples": samples, "quantified": quantified, "events": events, "complete": complete}


#
# Result tables
#


def suppa_event_type(event_id):
    match = SUPPA_EVENT_ID.match(event_id)
    return SUPPA_EVENT_TYPES[match.group(1)] if match else None


def suppa_diffsplice(path, alpha, dpsi_cutoff, genes):
    """
    Events of a SUPPA2 dPSI table, those tested (with a dPSI and a p-value) and those
    significant, overall, by direction and by event type.
    """
    calls = Calls()
    stats = {"events": 0, "significant_up": 0, "significant_down": 0, "tested_by_type": defaultdict(int), "significant_by_type": defaultdict(int)}
    with open(path) as handle:
        handle.readline()
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 3:
                continue
            stats["events"] += 1
            dpsi, pvalue = to_float(fields[1]), to_float(fields[2])
            if dpsi is None or pvalue is None:
                continue
            significant = pvalue < alpha and abs(dpsi) > dpsi_cutoff
            calls.add(significant, genes.key(fields[0].split(";")[0]), pvalue)
            event_type = suppa_event_type(fields[0])
            if event_type:
                stats["tested_by_type"][event_type] += 1
            if significant:
                stats["significant_up" if dpsi > 0 else "significant_down"] += 1
                if event_type:
                    stats["significant_by_type"][event_type] += 1
    stats["calls"] = calls
    return stats


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
        gene = genes.key(row["gene_id"])
        isoforms[contrast].add(significant)
        gene_switches[contrast][gene] |= significant
    gene_calls = {}
    for contrast, switches in gene_switches.items():
        gene_calls[contrast] = Calls()
        for gene, switched in switches.items():
            gene_calls[contrast].add(switched, gene)
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
                gene = genes.key(gene)
                if gene:
                    calls.tested_genes.add(gene)
                    if significant:
                        calls.significant_genes.add(gene)
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
        if fdr is None:
            continue
        calls.add(fdr <= RMATS_FDR, genes.key(row["GeneID"]), to_float(row.get("PValue")))


#
# Sections
#


def write_general_stats(general_stats):
    columns = [column for column in GENERAL_STATS_COLUMNS if any(column in values for values in general_stats.values())]
    with open(GENERAL_STATS_FILE, "w") as handle:
        handle.write("\\t".join(["Sample"] + columns) + "\\n")
        for sample in sorted(general_stats):
            values = general_stats[sample]
            cells = ["" if values.get(column) is None else GENERAL_STATS_COLUMNS[column].format(values[column]) for column in columns]
            handle.write("\\t".join([sample] + cells) + "\\n")


def write_dexseq_section(assignment):
    keys = ["in_bins"] + [key for key, _, _ in DEXSEQ_COUNTERS.values()]
    categories = {"in_bins": {"name": "Counted in exon bins", "color": "#437bb1"}}
    for key, name, colour in DEXSEQ_COUNTERS.values():
        if any(key in counts for counts in assignment.values()):
            categories[key] = {"name": name, "color": colour}
    write_section(
        "dexseq_read_assignment",
        "DEXSeq: Read Assignment",
        "Reads that `dexseq_count.py` counted into exon bins and the reads it could not assign, in each sample. "
        "The numbers are approximate: a read overlapping several exon bins of one gene is counted once for each bin.",
        "bargraph",
        {sample: {key: counts.get(key, 0) for key in keys if key in categories} for sample, counts in assignment.items()},
        helptext="""
            `dexseq_count.py` counts the reads overlapping each exon bin of the flattened
            annotation, and writes counters of the reads it could not assign:

            * `_empty`: the read overlaps no exon bin
            * `_ambiguous`: the read overlaps exon bins of more than one gene
            * `_lowaqual`: the alignment quality is below the minimum
            * `_notaligned`: the read is not aligned

            A read overlapping several exon bins of one gene is counted once for each bin, and
            the script does not write how many reads it counted, so the shares are approximate.
            A large `_empty` share often means `dexseq_count` ran with the wrong `--stranded`
            setting for the library.
            """,
        pconfig={
            "id": "dexseq_read_assignment_plot",
            "title": "DEXSeq: Read Assignment",
            "ylab": "Reads",
            "cpswitch_counts_label": "Number of reads",
        },
        categories=categories,
    )


def write_suppa_psi_sections(psi_by_route):
    route = first_route(psi_by_route)
    if route is None:
        return
    psi = psi_by_route[route]
    note = f" Both quantification routes ran, and this shows the {route_label(route)} route." if len(psi_by_route) > 1 else ""

    write_section(
        "suppa_psi_quantification",
        "SUPPA2: PSI Quantification",
        f"Number of local splicing events (`psiPerEvent`) for which SUPPA2 could compute a PSI value in each sample.{note}",
        "bargraph",
        {sample: {"quantified": psi["quantified"][index], "not_quantified": psi["events"] - psi["quantified"][index]} for index, sample in enumerate(psi["samples"])},
        helptext="""
            SUPPA2 reports `nan` instead of a PSI value when the transcripts that define an event
            are not expressed in a sample, when their total expression is below the
            `--psiperevent_total_filter` threshold, or when a transcript is missing from the
            expression file.

            A sample with a much lower fraction of quantified events than others usually has low
            depth or low transcript-level expression, and contributes fewer events to the
            differential splicing tests (`diffSplice` leaves out events with a missing value in
            more than the `--diffsplice_nan_threshold` fraction of the samples of a condition).
            """,
        pconfig={
            "id": "suppa_psi_quantification_plot",
            "title": "SUPPA2: PSI Quantification",
            "ylab": "Events",
            "cpswitch_counts_label": "Number of events",
        },
        categories={
            "quantified": {"name": "Quantified (PSI value)", "color": "#437bb1"},
            "not_quantified": {"name": "Not quantified (nan)", "color": "#b1b1b1"},
        },
    )

    # Events near 0 or 1 in every sample would dominate the correlation without saying anything about the samples
    informative = [values for values in psi["complete"] if 0.05 <= sum(values) / len(values) <= 0.95]
    if len(psi["samples"]) < 2 or len(informative) < 3:
        return
    step = max(1, math.ceil(len(informative) / PSI_CONCORDANCE_MAX_EVENTS))
    used = informative[::step]
    ranks = [centred_ranks([values[index] for values in used]) for index in range(len(psi["samples"]))]
    norms = [math.sqrt(sum(map(operator.mul, rank, rank))) for rank in ranks]
    matrix = [
        [round(sum(map(operator.mul, ranks[i], ranks[j])) / (norms[i] * norms[j]), 3) if norms[i] and norms[j] else None for j in range(len(ranks))]
        for i in range(len(ranks))
    ]
    sampled = f", {len(used):,} of them taken evenly across the matrix" if step > 1 else ""
    write_section(
        "suppa_psi_concordance",
        "SUPPA2: PSI Concordance",
        f"Spearman correlation of the PSI values of each pair of samples, over the {len(informative):,} local events "
        f"quantified in every sample with a mean PSI between 0.05 and 0.95{sampled}.{note}",
        "heatmap",
        matrix,
        helptext="""
            Samples of the same condition should correlate more closely with each other than with
            the samples of other conditions. A sample that correlates poorly with every other
            sample can point to a quantification or library problem, and one that clusters with
            the wrong condition to a sample swap.

            Events quantified in every sample are used, and those with a mean PSI below 0.05 or
            above 0.95 are left out, as they are nearly constant and would inflate every
            correlation.
            """,
        pconfig={
            "id": "suppa_psi_concordance_plot",
            "title": "SUPPA2: PSI Concordance",
            "xlab": "Sample",
            "ylab": "Sample",
            "max": 1,
            "square": True,
        },
        xcats=psi["samples"],
        ycats=psi["samples"],
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


def write_suppa_diffsplice_sections(diffsplice, significance):
    """diffsplice: (contrast, level, route) -> dPSI table summary."""
    if not diffsplice:
        return
    levels = {"events": "local events", "isoforms": "isoforms"}
    several_routes = len({route for _, _, route in diffsplice}) > 1
    keys = sorted(diffsplice, key=lambda key: (key[0], list(levels).index(key[1]), ROUTE_ORDER.index(key[2])))

    def comparison(contrast, level, route):
        name = contrast if level is None else f"{contrast} {levels[level]}"
        return f"{name} ({route_label(route)})" if several_routes and route else name

    rows = {}
    for contrast, level, route in keys:
        stats = diffsplice[(contrast, level, route)]
        calls = stats["calls"]
        row = {"events": stats["events"], "tested": calls.tested, "significant": calls.significant}
        row["significant_up"] = stats["significant_up"]
        row["significant_down"] = stats["significant_down"]
        if stats["events"]:
            row["pct_tested"] = 100 * calls.tested / stats["events"]
        if calls.tested:
            row["pct_significant"] = 100 * calls.significant / calls.tested
        rows[comparison(contrast, level, route)] = row

    count = {"format": "{:,.0f}"}
    direction = "dPSI is the PSI of condition 2 (the treatment) minus that of condition 1 (the control)"
    write_section(
        "suppa_differential_splicing",
        "SUPPA2: Differential Splicing",
        f"Local events and isoforms tested and called significant in each SUPPA2 `diffSplice` comparison. "
        f"Significant: {significance}. {direction}.",
        "table",
        rows,
        helptext=f"""
            Each row is one contrast of the contrastsheet, for the local events (`psiPerEvent`)
            or the isoforms (`psiPerIsoform`).

            Direction: the pipeline runs `diffSplice` with the control as condition 1 and the
            treatment as condition 2, and SUPPA2 reports condition 2 minus condition 1, so
            **Higher PSI** means a higher PSI in the treatment.

            Events are counted as significant with {significance}: the p-value threshold is
            `--diffsplice_alpha` and the dPSI threshold `--clusterevents_dpsithreshold`, which
            `clusterEvents` also uses to pick the events it clusters. `diffSplice` corrects its
            p-values for multiple testing by gene with `--diffsplice_gene_correction`.
            """,
        pconfig={
            "id": "suppa_differential_splicing_table",
            "title": "SUPPA2: Differential Splicing",
            "namespace": "SUPPA2",
            "col1_header": "Comparison",
            "sort_rows": False,
        },
        headers={
            "events": {**count, "title": "Events", "description": "Events in the dPSI table", "scale": "Greys", "hidden": True},
            "tested": {
                **count,
                "title": "Tested",
                "description": "Events with a dPSI value. Events with too many missing PSI values in a condition (see `--diffsplice_nan_threshold`) are not tested.",
                "scale": "Blues",
            },
            "pct_tested": {"title": "% Tested", "description": "Percentage of events with a dPSI value", "min": 0, "max": 100, "suffix": "%", "scale": "GnBu", "hidden": True},
            "significant": {**count, "title": "Significant", "description": f"Tested events with {significance}", "scale": "Purples"},
            "pct_significant": {
                "title": "% Significant",
                "description": f"Percentage of the tested events that are significant ({significance})",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdPu",
            },
            "significant_up": {
                **count,
                "title": "Higher PSI",
                "description": f"Significant events with dPSI > 0, a higher PSI in the treatment than in the control ({direction})",
                "scale": "Reds",
            },
            "significant_down": {
                **count,
                "title": "Lower PSI",
                "description": f"Significant events with dPSI < 0, a lower PSI in the treatment than in the control ({direction})",
                "scale": "Blues",
            },
        },
    )

    # Isoform comparisons have no event types, and the table already shows them
    local = {comparison(contrast, None, route): diffsplice[(contrast, level, route)] for contrast, level, route in keys if level == "events"}
    if not local:
        return
    with_significant = {name: stats for name, stats in local.items() if stats["calls"].significant}
    without_significant = [name for name, stats in local.items() if not stats["calls"].significant]
    description = (
        f"Significant local events (`psiPerEvent`) of each comparison by event type, as a percentage of the tested events "
        f"of the type, with the numbers of events in the second tab. Significant: {significance}."
    )
    helptext = """
        SUPPA2 event types:

        * **SE**: skipped exon
        * **A5** / **A3**: alternative 5' / 3' splice site
        * **MX**: mutually exclusive exons
        * **RI**: retained intron
        * **AF** / **AL**: alternative first / last exon

        Event types differ in how many events the annotation defines, so the percentage of the
        tested events of each type compares them more fairly than the numbers. Comparisons of
        isoforms (`psiPerIsoform`) are not shown, as they have no event types.
        """
    if not with_significant:
        write_section("suppa_significant_events", "SUPPA2: Significant Events", description, "html", "<p>No comparison has any significant events.</p>", helptext=helptext)
        return
    if without_significant:
        description += f" Comparisons without significant events, not shown: {', '.join(without_significant)}."
    categories = {event_type: {"name": event_type, "color": colour} for event_type, colour in SUPPA_EVENT_TYPE_COLOURS.items()}
    write_section(
        "suppa_significant_events",
        "SUPPA2: Significant Events",
        description,
        "bargraph",
        [
            {
                name: {
                    event_type: 100 * stats["significant_by_type"].get(event_type, 0) / tested
                    for event_type, tested in stats["tested_by_type"].items()
                    if tested
                }
                for name, stats in with_significant.items()
            },
            {name: dict(stats["significant_by_type"]) for name, stats in with_significant.items()},
        ],
        helptext=helptext,
        pconfig={
            "id": "suppa_significant_events_plot",
            "title": "SUPPA2: Significant Events",
            "ylab": "% of tested events of the type",
            "stacking": "group",
            "cpswitch": False,
            "data_labels": [
                {"name": "% of tested", "ylab": "% of tested events of the type", "ymax": 100},
                {"name": "Events", "ylab": "Significant events"},
            ],
        },
        categories=[categories, categories],
    )


def write_leafcutter_introns_section(introns):
    write_section(
        "leafcutter_introns_detected",
        "LeafCutter: Introns Detected",
        "Clustered introns with and without a supporting read in each sample.",
        "bargraph",
        {sample: {"detected": counts["detected"], "not_detected": counts["total"] - counts["detected"]} for sample, counts in introns.items()},
        helptext="""
            LeafCutter pools the junctions of all samples, keeps introns with enough reads across
            the cohort and groups introns that share a splice site into clusters. Every sample is
            counted against the same set of introns, and an intron is detected in a sample when at
            least one read in that sample supports it.

            A sample that detects a much smaller share of the introns than the others has lower
            depth or library complexity. `leafcutter_ds` tests these counts directly, so such a
            sample contributes little and can push clusters below its coverage filters. The
            junction reads of each sample are in the General Statistics table.
            """,
        pconfig={
            "id": "leafcutter_introns_detected_plot",
            "title": "LeafCutter: Introns Detected",
            "ylab": "Introns",
            "cpswitch_counts_label": "Number of introns",
        },
        categories={
            "detected": {"name": "Detected (at least one read)", "color": "#437bb1"},
            "not_detected": {"name": "Not detected", "color": "#b1b1b1"},
        },
    )


def write_leafcutter_outcomes_section(outcomes, args):
    categories = {
        "significant": {"name": f"Significant (adjusted p < {fmt(args.fdr)})", "color": "#7b4ea3"},
        "not_significant": {"name": "Not significant", "color": "#437bb1"},
        "no_pvalue": {"name": "Tested, no adjusted p-value", "color": "#a6c4e2"},
    }
    not_tested = [*LEAFCUTTER_SKIP_STATUSES.values(), ("timeout", "Not tested: timeout"), ("error", "Not tested: error")]
    for (key, label), colour in zip(not_tested, LEAFCUTTER_NOT_TESTED_COLOURS):
        categories[key] = {"name": label, "color": colour}
    settings = (
        f"`leafcutter_ds` ran with at least {args.leafcutter_min_samples_per_group} samples per group (`-g`, "
        f"`--leafcutter_ds_min_samples_per_group`), {args.leafcutter_min_samples_per_intron} samples per intron (`-i`, "
        f"`--leafcutter_ds_min_samples_per_intron`) and {args.leafcutter_min_coverage} reads per sample (`-c`, "
        f"`--leafcutter_ds_min_coverage`)."
    )
    write_section(
        "leafcutter_cluster_test_outcomes",
        "LeafCutter: Cluster Test Outcomes",
        f"Outcome of the differential splicing test for each intron cluster, one bar per contrast. "
        f"Significant: adjusted p < {fmt(args.fdr)}. {settings}",
        "bargraph",
        outcomes,
        helptext=f"""
            `leafcutter_ds` only tests clusters that pass its coverage filters. The reasons for
            not testing a cluster, as LeafCutter reports them in the `status` column, are:

            * **Too many introns**: `Too many introns in cluster` (more than `--max_cluster_size`)
            * **One junction**: `<=1 junction in cluster`
            * **One sample with reads**: `<=1 sample with coverage>0`
            * **One sample above min coverage**: `<=1 sample with coverage>min_coverage`
            * **Too few introns in enough samples**: `<2 introns used in >=min_samples_per_intron samples`
            * **Not enough valid samples**: fewer than `-g` samples in each group reach `-c` reads
            * **Timeout**: the model fit took longer than `--timeout`
            * **Error**: the model fit failed, and LeafCutter reported the error message

            A contrast with few tested clusters usually has too few samples per group or too
            little depth for the filters, rather than little differential splicing.

            The p-values are adjusted with the Benjamini-Hochberg method, and a tested cluster is
            significant with an adjusted p-value below {fmt(args.fdr)}.
            """,
        pconfig={
            "id": "leafcutter_cluster_test_outcomes_plot",
            "title": "LeafCutter: Cluster Test Outcomes",
            "ylab": "Clusters",
            "cpswitch_counts_label": "Number of clusters",
        },
        categories=categories,
    )


def reciprocal_contrasts(contrasts):
    by_pair = {(contrast["treatment"], contrast["control"]): name for name, contrast in contrasts.items()}
    return {name: by_pair[(contrast["control"], contrast["treatment"])] for name, contrast in contrasts.items() if (contrast["control"], contrast["treatment"]) in by_pair}


def write_differential_table(rows, reciprocal):
    """One row per contrast and one column per analysis, so the analyses of a contrast read
    across. Analyses count different features, so each column states its own in its header."""
    rows.sort(key=lambda row: (METHOD_ORDER.index(row["method"]), ROUTE_ORDER.index(row["route"]), row["unit"], row["contrast"]))
    headers = {}
    data = defaultdict(dict)
    if reciprocal:
        headers["reciprocal"] = {
            "title": "Reverse of",
            "description": "The contrast that compares the same conditions the other way round, with the same features called significant",
            "scale": False,
        }
        for contrast, other in reciprocal.items():
            data[contrast]["reciprocal"] = other
    for row in rows:
        key = re.sub("[^a-z0-9]+", "_", " ".join(part for part in (row["method"], row["route"], row["unit"]) if part).lower())
        route = f" ({route_label(row['route'])})" if row["route"] else ""
        title = f"{row['method']} {row['unit']}{route}"
        if key not in headers:
            route_text = f", from the {route_label(row['route'])} quantification" if row["route"] else ""
            headers[key] = {
                "title": title,
                "description": f"{row['method']}: {row['unit']} called significant ({row['threshold']}; {row['effect']}){route_text}{row['note']}",
                "format": "{:,.0f}",
                "scale": "Purples",
                "hidden": row["hidden"],
            }
            headers[f"{key}_tested"] = {
                "title": f"{title}: tested",
                "description": f"{row['method']}: {row['tested_text'] or row['unit'] + ' tested'}{route_text}",
                "format": "{:,.0f}",
                "scale": "Blues",
                "hidden": True,
            }
        data[row["contrast"]][key] = row["significant"]
        data[row["contrast"]][f"{key}_tested"] = row["tested"]
    reciprocal_text = (
        " Contrasts marked in the Reverse of column compare the same conditions the other way round, so their counts are the same calls, not additional ones."
        if reciprocal
        else ""
    )
    write_section(
        "differential_splicing",
        "Differential Splicing",
        "How many features each differential splicing analysis called significant in each contrast. Each column is one "
        "analysis, and its header gives what the analysis tests (events, exon bins, genes, clusters, isoforms or bins), "
        "the threshold and any effect-size cutoff, which are the pipeline parameters where the tool has them. Counts are "
        "not comparable across columns: the analyses count different units with different tests and effect-size "
        "cutoffs, so compare the genes they call in the method agreement sections instead. Hidden columns give the "
        "number of features each analysis tested. Analyses that depend on a quantification name it: STAR + Salmon for "
        f"Salmon quantification of the STAR alignments, Salmon for Salmon pseudo-alignment or Salmon results.{reciprocal_text}",
        "table",
        {contrast: data[contrast] for contrast in sorted(data)},
        pconfig={
            "id": "differential_splicing_table",
            "title": "Differential splicing",
            "namespace": "Differential splicing",
            "col1_header": "Contrast",
            "sort_rows": False,
        },
        headers=headers,
    )


def write_design_section(contrasts, conditions, ran, args):
    """Samples per condition of each contrast, against the minimum sample settings of the tools that ran."""
    if not contrasts or not conditions:
        return
    total = len({sample for samples in conditions.values() for sample in samples})
    headers = {
        "treatment": {"title": "Treatment", "description": "Treatment condition of the contrast", "scale": False},
        "control": {"title": "Control", "description": "Control condition of the contrast", "scale": False},
        "n_treatment": {"title": "Treatment samples", "description": "Samples of the treatment condition", "format": "{:,.0f}", "scale": False},
        "n_control": {"title": "Control samples", "description": "Samples of the control condition", "format": "{:,.0f}", "scale": False},
    }
    if "LeafCutter" in ran:
        headers["leafcutter_g"] = {
            "title": "LeafCutter min per group",
            "description": "Samples each group needs with enough reads for leafcutter_ds to test a cluster (-g, --leafcutter_ds_min_samples_per_group)",
            "format": "{:,.0f}",
            "scale": False,
        }
        headers["leafcutter_i"] = {
            "title": "LeafCutter min per intron",
            "description": "Samples an intron needs reads in for leafcutter_ds to use it (-i, --leafcutter_ds_min_samples_per_intron)",
            "format": "{:,.0f}",
            "scale": False,
        }
    if "DRIMSeq" in ran:
        headers["drimseq_gene"] = {
            "title": "DRIMSeq min samples",
            "description": "Samples a gene must be expressed in to pass the DRIMSeq filter, which runs once on all the samples of the run (--min_samps_gene_expr)",
            "format": "{:,.0f}",
            "scale": False,
        }
    headers["status"] = {
        "title": "Status",
        "description": "OK, or the minimum sample settings the contrast does not meet",
        "scale": False,
        "cond_formatting_rules": {"pass": [{"s_eq": "OK"}], "fail": [{"s_contains": "Below minimum"}], "warn": [{"s_contains": "Check"}]},
    }

    data = {}
    for name in sorted(contrasts):
        contrast = contrasts[name]
        n_treatment = len(conditions.get(contrast["treatment"], ()))
        n_control = len(conditions.get(contrast["control"], ()))
        row = {"treatment": contrast["treatment"], "control": contrast["control"], "n_treatment": n_treatment, "n_control": n_control}
        problems = []
        if min(n_treatment, n_control) < 2:
            problems.append("Below minimum: no replicates in a condition, so the tests have no within-condition variance")
        if "LeafCutter" in ran:
            row["leafcutter_g"] = args.leafcutter_min_samples_per_group
            row["leafcutter_i"] = args.leafcutter_min_samples_per_intron
            if min(n_treatment, n_control) < args.leafcutter_min_samples_per_group:
                problems.append(f"Below minimum: LeafCutter needs {args.leafcutter_min_samples_per_group} per group")
            if n_treatment + n_control < args.leafcutter_min_samples_per_intron:
                problems.append(f"Below minimum: LeafCutter needs {args.leafcutter_min_samples_per_intron} samples per intron")
        if "DRIMSeq" in ran:
            row["drimseq_gene"] = args.min_samps_gene_expr
            needed = max(args.min_samps_gene_expr, args.min_samps_feature_expr, args.min_samps_feature_prop)
            if total < needed:
                problems.append(f"Below minimum: DRIMSeq needs {needed} samples in the run")
        if "SUPPA2" in ran and args.suppa_paired and n_treatment != n_control:
            problems.append("Check: the SUPPA2 paired test needs the same number of samples in each condition")
        if "rMATS" in ran and args.rmats_paired_stats and n_treatment != n_control:
            problems.append("Check: rMATS paired statistics need the same number of samples in each condition")
        row["status"] = "; ".join(problems) or "OK"
        data[name] = row

    write_section(
        "splicing_design_check",
        "Design Check",
        "Samples of each condition of every contrast, next to the minimum sample settings of the analyses that ran. "
        "A contrast below a minimum gets few or no tested features from that analysis, however strong the splicing changes.",
        "table",
        data,
        pconfig={"id": "splicing_design_check_table", "title": "Splicing: Design Check", "namespace": "Design check", "col1_header": "Contrast", "sort_rows": False},
        headers=headers,
    )


def write_agreement_sections(gene_calls, genes):
    """gene_calls: contrast -> (method, route) -> Calls of genes."""
    usable = {contrast: methods for contrast, methods in gene_calls.items() if len({method for method, _ in methods}) >= 2}
    if usable:
        rows = {}
        capped = []
        for contrast in sorted(usable):
            families = defaultdict(set)
            for (method, _), calls in usable[contrast].items():
                for gene in calls.significant_genes:
                    families[gene].add(method)
            shared = sorted(((gene, methods) for gene, methods in families.items() if len(methods) >= 2), key=lambda item: (-len(item[1]), item[0]))
            if len(shared) > AGREEMENT_MAX_GENES:
                capped.append(f"{contrast} ({len(shared):,})")
            for gene, methods in shared[:AGREEMENT_MAX_GENES]:
                rows[f"{contrast}: {gene}"] = {
                    "contrast": contrast,
                    "gene_name": genes.name(gene) or gene,
                    "n_methods": len(methods),
                    "methods": ", ".join(method for method in METHOD_ORDER if method in methods),
                }
        description = (
            "Genes called significant by at least two of the analyses in each contrast. Each analysis calls a gene "
            "significant when it calls any of the gene's features significant, and the routes of one analysis count once."
        )
        if capped:
            description += f" Only the {AGREEMENT_MAX_GENES} genes with the most analyses are shown for: {', '.join(capped)}."
        if rows:
            write_section(
                "splicing_method_agreement",
                "Method Agreement: Shared Genes",
                description,
                "table",
                rows,
                helptext="""
                    The analyses test different features (exon bins, events, transcripts, isoforms, intron
                    clusters), so their calls are compared at the gene level. Gene IDs lose their version
                    suffix, and LeafCutter and IsoformSwitchAnalyzeR, which report gene names, are mapped to
                    gene IDs through the annotation where the name is unique.

                    A gene called by several analyses of different kinds is more likely a real change than one
                    called by a single analysis, but the analyses have different power: one that tests few
                    genes calls few, without disagreeing with the others.
                    """,
                pconfig={"id": "splicing_method_agreement_table", "title": "Splicing: Shared Genes", "namespace": "Method agreement", "col1_header": "Contrast: gene", "sort_rows": False},
                headers={
                    "contrast": {"title": "Contrast", "scale": False},
                    "gene_name": {"title": "Gene name", "scale": False},
                    "n_methods": {"title": "Analyses", "description": "Analyses that call the gene significant", "format": "{:,.0f}", "scale": "Purples"},
                    "methods": {"title": "Called by", "description": "Analyses that call the gene significant", "scale": False},
                },
            )
        else:
            write_section(
                "splicing_method_agreement",
                "Method Agreement: Shared Genes",
                description,
                "html",
                "<p>No gene is called significant by two or more analyses in any contrast.</p>",
            )

        labels_in_order = sorted(
            {key for methods in usable.values() for key in methods},
            key=lambda key: (METHOD_ORDER.index(key[0]), ROUTE_ORDER.index(key[1])),
        )
        xcats = [with_route(method, route) for method, route in labels_in_order]
        ycats = []
        matrix = []
        for contrast in sorted(usable):
            for key in labels_in_order:
                if key not in usable[contrast]:
                    continue
                ycats.append(f"{contrast}: {with_route(*key)}")
                matrix.append([jaccard(usable[contrast][key], usable[contrast].get(other)) for other in labels_in_order])
        write_section(
            "splicing_method_agreement_heatmap",
            "Method Agreement: Overlap",
            "Jaccard index of the genes each pair of analyses calls significant in each contrast, counting only the genes "
            "both analyses tested: the genes both call significant, over the genes either calls significant. Empty cells "
            "have no significant gene in either analysis.",
            "heatmap",
            matrix,
            pconfig={"id": "splicing_method_agreement_heatmap_plot", "title": "Splicing: Method Overlap", "xlab": "Analysis", "ylab": "Contrast: analysis", "min": 0, "max": 1},
            xcats=xcats,
            ycats=ycats,
        )

    routes = {}
    for contrast in sorted(gene_calls):
        methods = gene_calls[contrast]
        for method in METHOD_ORDER:
            if (method, "STAR") in methods and (method, "Salmon") in methods:
                star, salmon = methods[(method, "STAR")], methods[(method, "Salmon")]
                common = star.tested_genes & salmon.tested_genes
                star_significant = star.significant_genes & common
                salmon_significant = salmon.significant_genes & common
                routes[f"{contrast}: {method}"] = {
                    "tested": len(common),
                    "star": len(star_significant),
                    "salmon": len(salmon_significant),
                    "both": len(star_significant & salmon_significant),
                    "jaccard": jaccard(star, salmon),
                }
    if routes:
        count = {"format": "{:,.0f}", "scale": "Purples"}
        write_section(
            "splicing_route_agreement",
            "Method Agreement: Quantification Routes",
            "Genes the same analysis calls significant from the STAR + Salmon and the Salmon quantifications, over the genes "
            "both routes tested. The routes quantify the same samples, so most of their calls should agree.",
            "table",
            routes,
            pconfig={"id": "splicing_route_agreement_table", "title": "Splicing: Route Agreement", "namespace": "Route agreement", "col1_header": "Contrast: analysis", "sort_rows": False},
            headers={
                "tested": {**count, "title": "Tested by both", "description": "Genes both routes tested", "scale": "Blues"},
                "star": {**count, "title": "STAR + Salmon", "description": "Genes called significant from the STAR + Salmon quantification"},
                "salmon": {**count, "title": "Salmon", "description": "Genes called significant from the Salmon quantification"},
                "both": {**count, "title": "Both", "description": "Genes called significant from both quantifications"},
                "jaccard": {"title": "Jaccard", "description": "Genes called significant by both routes over genes called by either", "min": 0, "max": 1, "format": "{:.2f}", "scale": "RdYlGn"},
            },
        )


def jaccard(first, second):
    if second is None:
        return None
    common = first.tested_genes & second.tested_genes
    union = (first.significant_genes | second.significant_genes) & common
    if not union:
        return None
    return round(len(first.significant_genes & second.significant_genes & common) / len(union), 3)


def write_pvalue_section(histograms):
    """histograms: contrast -> series name -> Calls."""
    datasets = {contrast: {name: calls for name, calls in series.items() if calls.p_count} for contrast, series in histograms.items()}
    datasets = {contrast: series for contrast, series in datasets.items() if series}
    if not datasets:
        return
    contrasts = sorted(datasets)
    write_section(
        "splicing_pvalue_histograms",
        "P-value Distributions",
        f"Distribution of the p-values of the tested features of each analysis, in {P_BINS} bins, as a percentage of the "
        "features with a p-value. One tab per contrast.",
        "linegraph",
        [
            {
                name: {f"{(position + 0.5) / P_BINS:.3f}": 100 * count / calls.p_count for position, count in enumerate(calls.p_bins)}
                for name, calls in datasets[contrast].items()
            }
            for contrast in contrasts
        ],
        helptext=f"""
            Under the null hypothesis p-values are uniform, so each bin holds about {100 / P_BINS:g}% of
            them. A well-behaved test shows a flat distribution with a peak near 0 from the real
            changes. A peak near 1, or a hump in the middle, points to a test that is
            conservative or a model that fits the data poorly, and a distribution that rises
            towards 0 everywhere to an anti-conservative test or an unmodelled batch effect.

            IsoformSwitchAnalyzeR and stageR report only adjusted values, so they have no curve.
            """,
        pconfig={
            "id": "splicing_pvalue_histograms_plot",
            "title": "Splicing: P-value Distributions",
            "xlab": "p-value",
            "ylab": "% of tested features",
            "xmin": 0,
            "xmax": 1,
            "ymin": 0,
            "tt_suffix": "%",
            "data_labels": [{"name": contrast, "ylab": "% of tested features"} for contrast in contrasts],
        },
    )


def main():
    args = parse_ext_args("${task.ext.args ?: ''}")
    fdr = fmt(args.fdr)
    correction = "gene-corrected" if args.suppa_gene_correction else "uncorrected"
    suppa_significance = f"{correction} p < {fmt(args.suppa_alpha)} and |dPSI| > {fmt(args.suppa_dpsi)}"

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
    conditions = defaultdict(set)
    for _, path in by_kind.pop("samplesheet", []):
        for row in read_delimited(path, ","):
            conditions[row["condition"]].add(row["sample"])

    general_stats = defaultdict(dict)
    dexseq_assignment = {}
    suppa_psi_by_route = {}
    suppa_diffsplice_by_key = {}
    leafcutter_intron_counts = {}
    leafcutter_outcomes_by_contrast = {}
    gene_calls = defaultdict(dict)
    histograms = defaultdict(dict)
    rmats_calls = defaultdict(Calls)
    rows = []
    ran = set()

    def add_row(method, route, contrast, unit, counts, threshold, effect=NO_EFFECT_CUTOFF, tested_text=None, hidden=False, note=""):
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
                "hidden": hidden,
                "note": note,
            }
        )

    dexseq_tested = "after independent filtering, the {} with an adjusted p-value"

    for kind, entries in by_kind.items():
        for label, path in entries:
            # The id is the sample of a per-sample table and the contrast of a result table
            route, contrast = label["route"], label["id"]
            if kind == "star_log":
                total, annotated = star_splices(path)
                if total:
                    general_stats[contrast]["star_annotated_splices"] = percentage(annotated or 0, total)
                    general_stats[contrast]["star_splices"] = total / 1e6
            elif kind == "dexseq_counts":
                dexseq_assignment[contrast] = dexseq_read_assignment(path)
            elif kind == "suppa_psi":
                suppa_psi_by_route[route] = suppa_psi(path)
            elif kind == "leafcutter_numers":
                leafcutter_intron_counts.update(leafcutter_introns(path))
            elif kind == "dexseq_deu_bins":
                calls = scan(path, "padj", args.fdr, genes, p_column="pvalue")
                histograms[contrast]["DEXSeq DEU exon bins"] = calls
                add_row("DEXSeq DEU", route, contrast, "exon bins", calls.counts(), f"padj < {fdr}", tested_text=dexseq_tested.format("exon bins"))
            elif kind == "dexseq_deu_genes":
                calls = scan(path, "padj", args.fdr, genes, gene_column="groupID")
                gene_calls[contrast][("DEXSeq DEU", route)] = calls
                add_row("DEXSeq DEU", route, contrast, "genes", calls.counts(), f"gene q-value < {fdr}", tested_text=dexseq_tested.format("genes"))
            elif kind == "edger_exons":
                calls = scan(path, "FDR", args.fdr, genes, p_column="P.Value", delimiter=",")
                histograms[contrast]["edgeR DEU exons"] = calls
                add_row("edgeR DEU", route, contrast, "exons", calls.counts(), f"FDR < {fdr}")
            elif kind == "edger_genes":
                calls = scan(path, "FDR", args.fdr, genes, gene_column="Geneid", delimiter=",")
                gene_calls[contrast][("edgeR DEU", route)] = calls
                add_row("edgeR DEU", route, contrast, "genes", calls.counts(), f"Simes FDR < {fdr}")
            elif kind == "dexseq_dtu_transcripts":
                ran.add("DRIMSeq")
                calls = scan(path, "padj", args.fdr, genes, p_column="pvalue")
                histograms[contrast][with_route("DEXSeq DTU transcripts", route)] = calls
                add_row(
                    "DEXSeq DTU",
                    route,
                    contrast,
                    "transcripts",
                    calls.counts(),
                    f"padj < {fdr}",
                    tested_text=dexseq_tested.format("transcripts"),
                    hidden=True,
                    note=". Without stageR, transcript-level padj does not control the error rate across genes; the stageR column is the transcript call",
                )
            elif kind == "dexseq_dtu_genes":
                ran.add("DRIMSeq")
                calls = scan(path, "padj", args.fdr, genes, gene_column="groupID")
                add_row("DEXSeq DTU", route, contrast, "genes", calls.counts(), f"gene q-value < {fdr}", tested_text=dexseq_tested.format("genes"))
            elif kind == "stager":
                ran.add("DRIMSeq")
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
                ran.add("SUPPA2")
                level = kind.removeprefix("suppa_")
                stats = suppa_diffsplice(path, args.suppa_alpha, args.suppa_dpsi, genes)
                suppa_diffsplice_by_key[(contrast, level, route)] = stats
                unit = "events" if level == "events" else "transcripts"
                histograms[contrast][with_route(f"SUPPA2 {unit}", route)] = stats["calls"]
                if level == "events":
                    gene_calls[contrast][("SUPPA2", route)] = stats["calls"]
                add_row("SUPPA2", route, contrast, unit, stats["calls"].counts(), suppa_significance.rsplit(" and ", 1)[0], effect=f"|dPSI| > {fmt(args.suppa_dpsi)}")
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
                    )
            elif kind == "leafcutter_ds":
                ran.add("LeafCutter")
                outcomes, calls = leafcutter_outcomes(path, args.fdr, genes)
                leafcutter_outcomes_by_contrast[contrast] = outcomes
                tested = sum(outcomes.get(key, 0) for key in ("significant", "not_significant", "no_pvalue"))
                histograms[contrast]["LeafCutter clusters"] = calls
                gene_calls[contrast][("LeafCutter", route)] = calls
                add_row("LeafCutter", route, contrast, "clusters", (tested, outcomes.get("significant", 0)), f"p.adjust < {fdr}")
            elif kind == "aspli_bins":
                calls = scan(path, "bin.fdr", args.fdr, genes, gene_column="locus", p_column="pvalue")
                histograms[contrast]["ASpli bins"] = calls
                gene_calls[contrast][("ASpli", route)] = calls
                add_row("ASpli", route, contrast, "bins", calls.counts(), f"bin.fdr < {fdr}")
            elif kind == "aspli_clusters":
                calls = scan(path, "FDR", args.fdr, genes, p_column="pvalue")
                histograms[contrast]["ASpli junction clusters"] = calls
                add_row("ASpli", route, contrast, "junction clusters", calls.counts(), f"FDR < {fdr}")
            elif kind == "rmats_summary":
                ran.add("rMATS")
                add_row(
                    "rMATS",
                    route,
                    contrast,
                    "JC events",
                    rmats_summary(path),
                    f"FDR <= {fmt(RMATS_FDR)}",
                    effect="no |dPSI| cutoff, the rMATS default",
                    tested_text="events with junction-count (JC) statistics, the primary rMATS result",
                )
            elif kind == "rmats_jc":
                ran.add("rMATS")
                rmats_jc(path, rmats_calls[contrast], genes)
            else:
                raise ValueError(f"Unknown splicing result kind '{kind}' for {path}")

    for contrast, calls in rmats_calls.items():
        gene_calls[contrast][("rMATS", "")] = calls
        histograms[contrast]["rMATS JC events"] = calls

    for sample, counts in dexseq_assignment.items():
        general_stats[sample]["dexseq_reads_in_bins"] = percentage(counts["in_bins"], sum(counts.values()))
    for sample, counts in leafcutter_intron_counts.items():
        general_stats[sample]["leafcutter_junction_reads"] = counts["reads"]
        general_stats[sample]["leafcutter_introns_detected"] = percentage(counts["detected"], counts["total"])

    if general_stats:
        write_general_stats(general_stats)
    if rows:
        write_differential_table(rows, reciprocal_contrasts({name: contrast for name, contrast in contrasts.items() if any(row["contrast"] == name for row in rows)}))
    write_design_section(contrasts, conditions, ran, args)
    write_agreement_sections(gene_calls, genes)
    write_pvalue_section(histograms)
    if dexseq_assignment:
        write_dexseq_section(dexseq_assignment)
    write_suppa_psi_sections(suppa_psi_by_route)
    write_suppa_diffsplice_sections(suppa_diffsplice_by_key, suppa_significance)
    if leafcutter_intron_counts:
        write_leafcutter_introns_section(leafcutter_intron_counts)
    if leafcutter_outcomes_by_contrast:
        write_leafcutter_outcomes_section(leafcutter_outcomes_by_contrast, args)

    with open("versions.yml", "w") as handle:
        handle.write('"${task.process}":\\n')
        handle.write(f"    python: {platform.python_version()}\\n")


if __name__ == "__main__":
    main()
