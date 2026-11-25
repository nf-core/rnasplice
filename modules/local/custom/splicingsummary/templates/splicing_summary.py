#!/usr/bin/env python

# Summarise the per-event, per-intron, per-exon-bin and per-transcript tables of the
# splicing analyses into two small MultiQC custom-content tables: per-sample
# detection rates for the General Statistics table, and per-contrast counts of
# tested and significant features.

import argparse
import csv
import gzip
import json
import math
import platform
import shlex
from collections import defaultdict

labels = json.loads(r"""${groovy.json.JsonOutput.toJson(labels)}""")
staged = json.loads(r"""${groovy.json.JsonOutput.toJson((results instanceof List ? results : [results]).collect { result -> result.toString() })}""")
GENERAL_STATS_FILE = "splicing_sample_metrics_mqc.tsv"
TABLE_FILE = "differential_splicing_mqc.json"

GENERAL_STATS_COLUMNS = {
    "dexseq_counts": "dexseq_exon_bins_detected",
    "suppa_psi": "suppa_events_quantified",
    "leafcutter_numers": "leafcutter_introns_detected",
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
]

ROUTE_ORDER = ["STAR", "Salmon", ""]


def parse_ext_args(args_string):
    if args_string in (None, "", "null"):
        args_string = ""
    parser = argparse.ArgumentParser()
    parser.add_argument("--fdr", type=float, default=0.05, help="Adjusted p-value cutoff of the tools without one of their own")
    parser.add_argument("--suppa_alpha", type=float, default=0.05, help="p-value cutoff of SUPPA2 diffSplice")
    parser.add_argument("--isa_alpha", type=float, default=0.05, help="Switch q-value cutoff of IsoformSwitchAnalyzeR")
    parser.add_argument("--isa_dif", type=float, default=0.1, help="Absolute dIF cutoff of IsoformSwitchAnalyzeR")
    return parser.parse_args(shlex.split(args_string))


def to_float(value):
    """A number, or None for the missing values R and Python write."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, newline="")


def read_delimited(path, delimiter="\\t"):
    """Rows of a table with a header, as dictionaries, one at a time."""
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def count_tested(path, column, cutoff, delimiter="\\t"):
    tested = significant = 0
    for row in read_delimited(path, delimiter):
        value = to_float(row.get(column))
        if value is None:
            continue
        tested += 1
        if value < cutoff:
            significant += 1
    return tested, significant


def dexseq_bins_detected(path):
    """Percentage of exon bins with a read, leaving out the HTSeq counters."""
    bins = detected = 0
    with open(path) as handle:
        for line in handle:
            feature, _, count = line.rstrip("\\n").partition("\\t")
            if not feature or feature.startswith("_"):
                continue
            bins += 1
            if int(count) > 0:
                detected += 1
    return 100 * detected / bins if bins else None


def suppa_events_quantified(path):
    """Percentage of the events of each sample with a PSI, from a SUPPA2 PSI matrix."""
    with open(path) as handle:
        samples = handle.readline().rstrip("\\n").split("\\t")
        quantified = [0] * len(samples)
        events = 0
        for line in handle:
            values = line.rstrip("\\n").split("\\t")[1:]
            events += 1
            for index, value in enumerate(values):
                if to_float(value) is not None:
                    quantified[index] += 1
    return {sample: 100 * quantified[index] / events for index, sample in enumerate(samples)} if events else {}


def leafcutter_introns_detected(path):
    """Percentage of the clustered introns of each sample with a read."""
    with open_text(path) as handle:
        samples = handle.readline().split()
        detected = [0] * len(samples)
        introns = 0
        for line in handle:
            counts = line.split()[1:]
            introns += 1
            for index, count in enumerate(counts):
                if int(count) > 0:
                    detected[index] += 1
    return {sample: 100 * detected[index] / introns for index, sample in enumerate(samples)} if introns else {}


def suppa_dpsi(path, alpha):
    """Events with both a dPSI and a p-value are tested."""
    tested = significant = 0
    with open(path) as handle:
        handle.readline()
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            dpsi, pvalue = to_float(fields[1]), to_float(fields[2])
            if dpsi is None or pvalue is None:
                continue
            tested += 1
            if pvalue < alpha:
                significant += 1
    return tested, significant


def leafcutter_clusters(path, fdr):
    tested = significant = 0
    for row in read_delimited(path):
        if row["status"] != "Success":
            continue
        tested += 1
        padj = to_float(row["p.adjust"])
        if padj is not None and padj < fdr:
            significant += 1
    return tested, significant


def stager_levels(path, fdr):
    """
    Genes from the screening stage, transcripts from the confirmation stage. stageR
    only confirms the transcripts of the genes that pass screening, so every
    transcript of a screened gene counts as tested.
    """
    genes = {}
    tested_tx = significant_tx = 0
    for row in read_delimited(path):
        gene_padj = to_float(row["gene"])
        if gene_padj is None:
            continue
        genes[row["geneID"]] = gene_padj
        tested_tx += 1
        tx_padj = to_float(row["transcript"])
        if tx_padj is not None and tx_padj < fdr:
            significant_tx += 1
    significant_genes = sum(1 for padj in genes.values() if padj < fdr)
    return (len(genes), significant_genes), (tested_tx, significant_tx)


def isoformswitchanalyzer(path, contrasts_path, alpha, dif):
    """Isoforms with a switch q-value are tested, by the contrast of their comparison."""
    contrasts = {}
    if contrasts_path:
        for row in read_delimited(contrasts_path, ","):
            contrasts[(row["control"], row["treatment"])] = row["contrast"]
    counts = defaultdict(lambda: [0, 0])
    for row in read_delimited(path, ","):
        qvalue = to_float(row["isoform_switch_q_value"])
        if qvalue is None:
            continue
        pair = (row["condition_1"], row["condition_2"])
        contrast = contrasts.get(pair, f"{pair[1]}-{pair[0]}")
        counts[contrast][0] += 1
        delta = to_float(row["dIF"])
        if qvalue < alpha and delta is not None and abs(delta) > dif:
            counts[contrast][1] += 1
    return counts


def write_differential_table(rows):
    """One row per contrast and one column per analysis, so the analyses of a contrast read
    across. Analyses count different features, so each column states its own in its header."""
    rows.sort(key=lambda row: (METHOD_ORDER.index(row["method"]), ROUTE_ORDER.index(row["route"]), row["unit"], row["contrast"]))
    headers = {}
    data = defaultdict(dict)
    for row in rows:
        key = "_".join(part for part in (row["method"], row["route"], row["unit"]) if part)
        key = "".join(char if char.isalnum() else "_" for char in key.lower())
        route = f" ({row['route']})" if row["route"] else ""
        title = f"{row['method']} {row['unit']}{route}"
        if key not in headers:
            route_text = f", from the {row['route']} quantification" if row["route"] else ""
            headers[key] = {
                "title": title,
                "description": f"{row['method']}: {row['unit']} called significant ({row['threshold']}){route_text}",
                "format": "{:,.0f}",
                "scale": "Purples",
            }
            headers[f"{key}_tested"] = {
                "title": f"{title}: tested",
                "description": f"{row['method']}: {row['unit']} tested{route_text}",
                "format": "{:,.0f}",
                "scale": "Blues",
                "hidden": True,
            }
        data[row["contrast"]][key] = row["significant"]
        data[row["contrast"]][f"{key}_tested"] = row["tested"]
    with open(TABLE_FILE, "w") as handle:
        json.dump({"id": "differential_splicing", "headers": headers, "data": {contrast: data[contrast] for contrast in sorted(data)}}, handle, indent=4)


def fmt(number):
    return f"{number:g}"


def main():
    args = parse_ext_args("${task.ext.args ?: ''}")
    fdr = fmt(args.fdr)

    general_stats = defaultdict(dict)
    suppa_psi_routes = {}
    rows = []
    contrasts_path = None

    def add_row(method, route, contrast, unit, counts, threshold):
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
            }
        )

    for label, path in zip(labels, staged):
        if label["kind"] == "contrasts":
            contrasts_path = path

    for label, path in zip(labels, staged):
        # The id is the sample of a per-sample table and the contrast of a result table
        kind, route, contrast = label["kind"], label["route"], label["id"]
        if kind == "dexseq_counts":
            general_stats[contrast]["dexseq_exon_bins_detected"] = dexseq_bins_detected(path)
        elif kind == "suppa_psi":
            suppa_psi_routes[route] = path
        elif kind == "leafcutter_numers":
            for sample, value in leafcutter_introns_detected(path).items():
                general_stats[sample]["leafcutter_introns_detected"] = value
        elif kind == "dexseq_deu_bins":
            add_row("DEXSeq DEU", route, contrast, "exon bins", count_tested(path, "padj", args.fdr), f"padj < {fdr}")
        elif kind == "dexseq_deu_genes":
            add_row("DEXSeq DEU", route, contrast, "genes", count_tested(path, "padj", args.fdr), f"gene q-value < {fdr}")
        elif kind == "edger_exons":
            add_row("edgeR DEU", route, contrast, "exons", count_tested(path, "FDR", args.fdr, ","), f"FDR < {fdr}")
        elif kind == "edger_genes":
            add_row("edgeR DEU", route, contrast, "genes", count_tested(path, "FDR", args.fdr, ","), f"Simes FDR < {fdr}")
        elif kind == "dexseq_dtu_transcripts":
            add_row("DEXSeq DTU", route, contrast, "transcripts", count_tested(path, "padj", args.fdr), f"padj < {fdr}")
        elif kind == "dexseq_dtu_genes":
            add_row("DEXSeq DTU", route, contrast, "genes", count_tested(path, "padj", args.fdr), f"gene q-value < {fdr}")
        elif kind == "stager":
            genes, transcripts = stager_levels(path, args.fdr)
            add_row("DEXSeq DTU + stageR", route, contrast, "genes", genes, f"screening padj < {fdr}")
            add_row("DEXSeq DTU + stageR", route, contrast, "transcripts", transcripts, f"confirmation padj < {fdr}")
        elif kind == "suppa_events":
            add_row("SUPPA2", route, contrast, "events", suppa_dpsi(path, args.suppa_alpha), f"p < {fmt(args.suppa_alpha)}")
        elif kind == "suppa_isoforms":
            add_row("SUPPA2", route, contrast, "transcripts", suppa_dpsi(path, args.suppa_alpha), f"p < {fmt(args.suppa_alpha)}")
        elif kind == "isa_features":
            threshold = f"q < {fmt(args.isa_alpha)} and |dIF| > {fmt(args.isa_dif)}"
            for isa_contrast, counts in isoformswitchanalyzer(path, contrasts_path, args.isa_alpha, args.isa_dif).items():
                add_row("IsoformSwitchAnalyzeR", route, isa_contrast, "isoforms", counts, threshold)
        elif kind == "leafcutter_ds":
            add_row("LeafCutter", route, contrast, "clusters", leafcutter_clusters(path, args.fdr), f"p.adjust < {fdr}")
        elif kind == "aspli_bins":
            add_row("ASpli", route, contrast, "bins", count_tested(path, "bin.fdr", args.fdr), f"bin.fdr < {fdr}")
        elif kind == "aspli_clusters":
            add_row("ASpli", route, contrast, "junction clusters", count_tested(path, "FDR", args.fdr), f"FDR < {fdr}")
        elif kind != "contrasts":
            raise ValueError(f"Unknown splicing result kind '{kind}' for {path}")

    # The STAR and Salmon routes quantify the same samples, so one matrix is enough
    for route in ROUTE_ORDER:
        if route in suppa_psi_routes:
            for sample, value in suppa_events_quantified(suppa_psi_routes[route]).items():
                general_stats[sample]["suppa_events_quantified"] = value
            break

    if general_stats:
        columns = [column for column in GENERAL_STATS_COLUMNS.values() if any(column in values for values in general_stats.values())]
        with open(GENERAL_STATS_FILE, "w") as handle:
            handle.write("\\t".join(["Sample"] + columns) + "\\n")
            for sample in sorted(general_stats):
                values = general_stats[sample]
                cells = ["" if values.get(column) is None else f"{values[column]:.2f}" for column in columns]
                handle.write("\\t".join([sample] + cells) + "\\n")

    if rows:
        write_differential_table(rows)

    with open("versions.yml", "w") as handle:
        handle.write('"${task.process}":\\n')
        handle.write(f"    python: {platform.python_version()}\\n")


if __name__ == "__main__":
    main()
