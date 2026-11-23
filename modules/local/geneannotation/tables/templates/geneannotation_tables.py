#!/usr/bin/env python3
"""
Write the gene and transcript annotation tables of a gffread transcript table.

The input is the table `gffread --table` writes with the columns transcript ID, gene ID,
gene name, `gene_biotype` and `gene_type`, in that order and with no header. gffread
writes `.` for an attribute a transcript does not have.

The gene table has one row per gene, with its name and biotype, and is the one the
results tables of the pipeline are annotated from: the results are by gene, and a
transcript table would give a gene one row per transcript in a join. The transcript
table adds the gene of each transcript to the same columns.

A gene with no name is named by its ID, so that every gene has a name to show. The
biotype is read from `gene_biotype` (Ensembl) or `gene_type` (GENCODE).

NOTE: This file is a Nextflow template. Nextflow interpolates the process variables
      and treats the backslash as an escape character, therefore every literal
      backslash must be doubled.
"""

import platform
import sys

MISSING = ["", "."]


def value(field):
    """Return a field of the table, or an empty string if gffread wrote it as missing."""
    field = field.strip()
    return "" if field in MISSING else field


def read_table(table):
    """Return the name and biotype of each gene, and the gene of each transcript."""
    genes = {}
    transcripts = {}
    with open(table) as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\\r\\n")
            if not line:
                continue
            fields = line.split("\\t")
            if len(fields) != 5:
                sys.exit(
                    f"Line {number} of {table} has {len(fields)} columns, where 5 are expected: "
                    "transcript ID, gene ID, gene name, gene_biotype and gene_type."
                )
            transcript_id, gene_id, gene_name, gene_biotype, gene_type = (value(field) for field in fields)
            if not gene_id:
                continue
            gene = genes.setdefault(gene_id, {"gene_name": "", "gene_biotype": ""})
            if not gene["gene_name"]:
                gene["gene_name"] = gene_name
            if not gene["gene_biotype"]:
                gene["gene_biotype"] = gene_biotype or gene_type
            if transcript_id:
                transcripts.setdefault(transcript_id, gene_id)
    return genes, transcripts


def main(table, prefix):
    genes, transcripts = read_table(table)
    if not genes:
        sys.exit(f"No transcript of {table} has a gene ID.")

    unnamed = [gene_id for gene_id, gene in genes.items() if not gene["gene_name"]]
    if unnamed:
        print(
            f"WARNING: {len(unnamed)} of {len(genes)} genes have no name in {table} and are named by their ID.",
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
    main("${tx2gene}", prefix)

    with open("versions.yml", "w") as handle:
        handle.write(f'"${task.process}":\\n    python: "{platform.python_version()}"\\n')
