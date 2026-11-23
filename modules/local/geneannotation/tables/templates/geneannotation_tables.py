#!/usr/bin/env python3
"""
Write the gene and transcript annotation tables of a GTF file.

The gene table has one row per gene, with its name and biotype, and is the one the
results tables of the pipeline are annotated from: the results are by gene, and a
transcript table would give a gene one row per transcript in a join. The transcript
table adds the gene of each transcript to the same columns.

A gene with no name in the GTF is named by its ID, so that every gene has a name to
show. The biotype is read from `gene_biotype` (Ensembl) or `gene_type` (GENCODE).

NOTE: This file is a Nextflow template. Nextflow interpolates the process variables
      and treats the backslash as an escape character, therefore every literal
      backslash must be doubled.
"""

import argparse
import gzip
import platform
import re
import shlex
import sys

ATTRIBUTE = re.compile(r'(\\S+)\\s+"([^"]*)"')
BIOTYPE_ATTRIBUTES = ["gene_biotype", "gene_type"]


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_id_attribute", default="gene_id")
    parser.add_argument("--gene_name_attribute", default="gene_name")
    return parser.parse_args(shlex.split(args))


def read_gtf(gtf, gene_id_attribute, gene_name_attribute):
    """Return the name and biotype of each gene, and the gene of each transcript."""
    genes = {}
    transcripts = {}
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 9:
                continue
            attributes = {}
            for key, value in ATTRIBUTE.findall(fields[8]):
                # The first value of an attribute, as a tag can be given more than once
                attributes.setdefault(key, value)
            gene_id = attributes.get(gene_id_attribute)
            if not gene_id:
                continue
            gene = genes.setdefault(gene_id, {"gene_name": "", "gene_biotype": ""})
            if not gene["gene_name"]:
                gene["gene_name"] = attributes.get(gene_name_attribute, "")
            if not gene["gene_biotype"]:
                gene["gene_biotype"] = next((attributes[key] for key in BIOTYPE_ATTRIBUTES if attributes.get(key)), "")
            transcript_id = attributes.get("transcript_id")
            if transcript_id:
                transcripts.setdefault(transcript_id, gene_id)
    return genes, transcripts


def main(gtf, prefix, args):
    options = parse_args(args)
    genes, transcripts = read_gtf(gtf, options.gene_id_attribute, options.gene_name_attribute)
    if not genes:
        sys.exit(f"No feature of {gtf} has a '{options.gene_id_attribute}' attribute.")

    unnamed = [gene_id for gene_id, gene in genes.items() if not gene["gene_name"]]
    if unnamed:
        print(
            f"WARNING: {len(unnamed)} of {len(genes)} genes have no '{options.gene_name_attribute}' in {gtf} and are named by their ID.",
            file=sys.stderr,
        )
    for gene_id in unnamed:
        genes[gene_id]["gene_name"] = gene_id
    for gene in genes.values():
        gene["gene_biotype"] = gene["gene_biotype"] or "NA"

    with open(f"{prefix}.gene_annotation.tsv", "w") as handle:
        handle.write("gene_id\\tgene_name\\tgene_biotype\\n")
        for gene_id in sorted(genes):
            gene = genes[gene_id]
            handle.write(f"{gene_id}\\t{gene['gene_name']}\\t{gene['gene_biotype']}\\n")

    with open(f"{prefix}.transcript_annotation.tsv", "w") as handle:
        handle.write("transcript_id\\tgene_id\\tgene_name\\tgene_biotype\\n")
        for transcript_id in sorted(transcripts):
            gene_id = transcripts[transcript_id]
            gene = genes[gene_id]
            handle.write(f"{transcript_id}\\t{gene_id}\\t{gene['gene_name']}\\t{gene['gene_biotype']}\\n")


if __name__ == "__main__":
    prefix = "${task.ext.prefix}" if "${task.ext.prefix}" not in ["null", ""] else "${meta.id}"
    main("${gtf}", prefix, "${task.ext.args ?: ''}")

    with open("versions.yml", "w") as handle:
        handle.write(f'"${task.process}":\\n    python: "{platform.python_version()}"\\n')
