#!/usr/bin/env python3
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import logging
import platform
from itertools import groupby
from typing import Iterator, Set

import yaml

logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("filter_gtf_for_genes_in_genome")
logger.setLevel(logging.INFO)


def is_header(line: str) -> bool:
    return line.startswith(">")


def extract_fasta_seq_names(fasta_name: str) -> Iterator[str]:
    """
    Yield the sequence names of a FASTA file.

    Modified from Brent Pedersen, "Correct Way To Parse A Fasta File In Python"
    https://www.biostars.org/p/710/
    """
    with open(fasta_name) as fh:
        # Group consecutive header/sequence lines together. Only the header
        # groups are of interest, the sequence groups are discarded.
        for is_header_group, group in groupby(fh, is_header):
            if not is_header_group:
                continue
            for line in group:
                # Drop the ">" and keep the first whitespace separated field
                fields = line[1:].strip().split()
                if not fields:
                    raise ValueError(f"Empty sequence name in FASTA file: {fasta_name}")
                yield fields[0]


def extract_genes_in_genome(fasta: str, gtf_in: str, gtf_out: str) -> None:
    seq_names_in_genome: Set[str] = set(extract_fasta_seq_names(fasta))
    logger.info("Extracted chromosome sequence names from : %s", fasta)
    logger.info("All chromosome names: %s", ", ".join(sorted(seq_names_in_genome)))

    seq_names_in_gtf: Set[str] = set()
    n_total_lines = 0
    n_lines_in_genome = 0

    with open(gtf_out, "w") as f, open(gtf_in) as g:
        for line in g:
            n_total_lines += 1
            seq_name_gtf = line.split("\\t")[0]
            seq_names_in_gtf.add(seq_name_gtf)
            if seq_name_gtf in seq_names_in_genome:
                n_lines_in_genome += 1
                f.write(line)

    logger.info(
        "Extracted %d / %d lines from %s matching sequences in %s",
        n_lines_in_genome,
        n_total_lines,
        gtf_in,
        fasta,
    )
    logger.info("All sequence IDs from GTF: %s", ", ".join(sorted(seq_names_in_gtf)))

    if n_lines_in_genome == 0:
        raise ValueError(
            f"No GTF lines in {gtf_in} match a sequence name in {fasta}. "
            "Please check that the GTF and the genome FASTA use the same naming convention."
        )

    logger.info("Wrote matching lines to %s", gtf_out)


if __name__ == "__main__":
    prefix = "${task.ext.prefix}" if "${task.ext.prefix}" not in ["null", ""] else "${fasta.baseName}"

    extract_genes_in_genome("${fasta}", "${gtf}", f"{prefix}_genes.gtf")

    # Write version information to versions.yml
    versions = {"${task.process}": {"python": platform.python_version()}}

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
