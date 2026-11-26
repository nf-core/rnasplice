#!/usr/bin/env python3
"""
Write a copy of a results table with the gene annotation next to its gene IDs.

The annotation columns (the gene name and biotype) are put right after the ID column,
and every row of the table is kept: a gene that is not in the annotation gets NA. The
gene ID of a row is looked up as it is, and if it is not found, without the version
suffix (`ENSG00000004961.15`) on either side, so that a table and an annotation that
differ only in their versions still match.

Two kinds of ID hold more than a gene ID:

- DEXSeq aggregates overlapping genes into one group, `ENSG00000147140+ENSG00000255591`,
  and SUPPA pools them as `ENSG00000147140_and_ENSG00000255591`. Each gene is looked up
  on its own, and the names are joined in the same way, `NONO+AL590762.1`.
- The SUPPA event IDs start with the gene ID, `ENSG00000004961;A5:X:273158-...`. Every
  ID is cut at its first `;` before it is looked up, as a GTF gene ID cannot hold one.

A header with one field fewer than the rows, as SUPPA and R write when the first column
holds row names, gets the ID column added to it.

The copy of a CSV table is a CSV, and of any other table a TSV.

NOTE: This file is a Nextflow template. Nextflow interpolates the process variables
      and treats the backslash as an escape character, therefore every literal
      backslash must be doubled.
"""

import csv
import os
import platform
import re
import sys

AGGREGATE_SEPARATOR = re.compile(r"(\\+|_and_)")
VERSION_SUFFIX = re.compile(r"\\.\\d+\\Z")


class GeneAnnotation:
    """The annotation of each gene, looked up with or without its version suffix."""

    def __init__(self, path):
        with open(path) as handle:
            rows = [line.rstrip("\\n").split("\\t") for line in handle if line.strip()]
        self.columns = rows[0][1:]
        self.genes = {row[0]: row[1:] for row in rows[1:]}
        # IDs without their version, for the IDs that lose it to only one gene
        unversioned = {}
        for gene_id in self.genes:
            unversioned.setdefault(VERSION_SUFFIX.sub("", gene_id), []).append(gene_id)
        self.unversioned = {key: ids[0] for key, ids in unversioned.items() if len(ids) == 1}

    def gene(self, gene_id):
        if gene_id in self.genes:
            return self.genes[gene_id]
        match = self.unversioned.get(VERSION_SUFFIX.sub("", gene_id))
        return self.genes[match] if match else None

    def annotate(self, gene_id):
        """The annotation of a gene ID, or of each gene of an aggregate, or None."""
        gene = self.gene(gene_id)
        if gene is not None:
            return gene
        parts = AGGREGATE_SEPARATOR.split(gene_id)
        if len(parts) == 1:
            return None
        genes = [self.gene(part) for part in parts[::2]]
        if all(gene is None for gene in genes):
            return None
        separators = parts[1::2] + [""]
        return [
            "".join((gene[column] if gene is not None else "NA") + separator for gene, separator in zip(genes, separators))
            for column in range(len(self.columns))
        ]


def main(table, annotation_path, id_column, prefix):
    annotation = GeneAnnotation(annotation_path)

    is_csv = table.endswith(".csv")
    output = f"{prefix}.annotated.{'csv' if is_csv else 'tsv'}"

    # A CSV can quote its fields, a TSV is split on its tabs as it is
    with open(table, newline="") as handle:
        if is_csv:
            rows = list(csv.reader(handle))
        else:
            rows = [line.rstrip("\\r\\n").split("\\t") for line in handle if line.strip()]

    # An empty table, as a stub or a tool can write, is copied as it is
    if not rows:
        open(output, "w").close()
        return

    header = rows[0]
    if len(rows) > 1 and len(header) == len(rows[1]) - 1:
        header = [id_column] + header
    if id_column not in header:
        sys.exit(f"{table} has no '{id_column}' column, its columns are: {', '.join(header)}")
    clashes = [column for column in annotation.columns if column in header]
    if clashes:
        sys.exit(f"{table} already has the column(s) {', '.join(clashes)} of the annotation.")
    id_index = header.index(id_column)

    unmatched = 0
    annotated = [header[: id_index + 1] + annotation.columns + header[id_index + 1 :]]
    for row in rows[1:]:
        gene_id = row[id_index].split(";", 1)[0]
        gene = annotation.annotate(gene_id)
        if gene is None:
            unmatched += 1
            gene = ["NA"] * len(annotation.columns)
        annotated.append(row[: id_index + 1] + gene + row[id_index + 1 :])

    with open(output, "w", newline="") as handle:
        if is_csv:
            csv.writer(handle, lineterminator="\\n").writerows(annotated)
        else:
            handle.writelines("\\t".join(row) + "\\n" for row in annotated)

    if unmatched:
        print(
            f"WARNING: {unmatched} of {len(rows) - 1} rows of {table} have a gene ID that is not in the annotation, and are annotated with NA.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    table = "${table}"
    prefix = "${task.ext.prefix}" if "${task.ext.prefix}" not in ["null", ""] else os.path.splitext(table)[0]
    main(table, "${gene_annotation}", "${id_column}", prefix)

    with open("versions.yml", "w") as handle:
        handle.write(f'"${task.process}":\\n    python: "{platform.python_version()}"\\n')
