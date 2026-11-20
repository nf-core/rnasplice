#!/usr/bin/env python3
"""
Give a strand to the junctions regtools could not strand.

An unstranded library has no read orientation to take the junction strand from, and
regtools writes `?` in the strand column unless the aligner left an `XS` tag behind.
The strand of a junction is however written in the genome itself: an intron whose
first and last two bases spell one of the canonical splice motifs can only be read in
one direction. Junctions with a non-canonical motif are looked up in the annotation,
and the ones that are in neither are left as they are, for the clustering to discard.

How many junctions, and how many reads, each of these routes gave a strand to is written
to a small table for MultiQC, so that what the clustering drops can be seen per sample.
"""

import os
import platform
import re
import yaml

# Intron boundary dinucleotides, and the strand each pair can only be read on. The
# reverse complement of a plus strand motif is the same motif seen on the minus strand
CANONICAL_MOTIFS = {
    ("GT", "AG"): "+",
    ("CT", "AC"): "-",
    ("GC", "AG"): "+",
    ("CT", "GC"): "-",
    ("AT", "AC"): "+",
    ("GT", "AT"): "-",
}

TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')

# Where the strand of each junction came from, and how the summary table names it
ORIGINS = {
    "already stranded": "already stranded (XS)",
    "motif": "stranded by motif",
    "annotation": "stranded by annotation",
    "strandless": "strandless",
}


class IndexedFasta:
    """Read pieces of a FASTA file through its `.fai` index, without loading it."""

    def __init__(self, fasta, fai):
        self.handle = open(fasta, "rb")
        self.index = {}
        with open(fai) as handle:
            for line in handle:
                name, length, offset, linebases, linewidth = line.split("\\t")[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))

    def fetch(self, chrom, start, end):
        """Return the sequence of `chrom` from `start` to `end`, both 1 based and included."""
        if chrom not in self.index:
            return ""
        length, offset, linebases, linewidth = self.index[chrom]
        if start < 1 or end > length:
            return ""
        # The line breaks are part of the file, so both the offset of the first base and
        # the number of bytes to read have to count them in
        first = offset + (start - 1) // linebases * linewidth + (start - 1) % linebases
        last = offset + (end - 1) // linebases * linewidth + (end - 1) % linebases
        self.handle.seek(first)
        return re.sub(r"[^A-Za-z]", "", self.handle.read(last - first + 1).decode()).upper()


def annotated_introns(gtf):
    """Map every intron between two consecutive exons of a transcript to its strand."""
    exons = {}
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            transcript = TRANSCRIPT_ID.search(fields[8])
            if not transcript:
                continue
            exons.setdefault(transcript.group(1), []).append(
                (fields[0], int(fields[3]), int(fields[4]), fields[6])
            )

    introns = {}
    for transcript_exons in exons.values():
        transcript_exons.sort(key=lambda exon: exon[1])
        for left, right in zip(transcript_exons, transcript_exons[1:]):
            chrom, _, left_end, strand = left
            intron = (chrom, left_end + 1, right[1] - 1)
            if intron[1] > intron[2]:
                continue
            # A transcript on either strand is no evidence of a strand at all
            if introns.get(intron, strand) != strand:
                introns[intron] = None
            else:
                introns[intron] = strand
    return introns


def intron_bounds(fields):
    """Take the intron a junction record spans, 1 based and both ends included.

    The record is the BED12 regtools writes: the two blocks are the anchors on each
    side of the intron, so the intron runs from the end of the first to the start of
    the second.
    """
    sizes = [int(size) for size in fields[10].rstrip(",").split(",")]
    if len(sizes) != 2:
        return None
    return (int(fields[1]) + sizes[0] + 1, int(fields[2]) - sizes[1])


def strand_junction(fields, genome, introns):
    """Give a strand to one junction record, and return where the strand came from."""
    if fields[5] != "?":
        return "already stranded"

    bounds = intron_bounds(fields)
    if bounds is None:
        return "strandless"

    start, end = bounds
    motif = (genome.fetch(fields[0], start, start + 1), genome.fetch(fields[0], end - 1, end))
    strand = CANONICAL_MOTIFS.get(motif)
    if strand:
        origin = "motif"
    else:
        strand = introns.get((fields[0], start, end))
        origin = "annotation"
    if not strand:
        return "strandless"

    fields[5] = strand
    return origin


def percent(part, total):
    return f"{100 * part / total:.2f}" if total else "NA"


def write_summary(sample, junctions, reads, output):
    """Write the number of junctions and of reads per origin of their strand, for MultiQC.

    The reads are the score column of the regtools records, which is the number of reads
    that span the junction, so a junction seen in many reads weighs more than a rare one.
    """
    header = ["Sample"]
    row = [sample]
    for unit, counts in (("Junctions", junctions), ("Reads", reads)):
        total = sum(counts.values())
        header += [unit] + [f"{unit} {label}" for label in ORIGINS.values()] + [f"{unit} strandless (%)"]
        row += [total] + [counts[origin] for origin in ORIGINS] + [percent(counts["strandless"], total)]

    with open(output, "w") as handle:
        handle.write("\\t".join(header) + "\\n")
        handle.write("\\t".join(str(value) for value in row) + "\\n")


def main(sample, junc, fasta, fai, gtf, output, summary):
    genome = IndexedFasta(fasta, fai)
    introns = annotated_introns(gtf)
    junctions = dict.fromkeys(ORIGINS, 0)
    reads = dict.fromkeys(ORIGINS, 0)

    with open(junc) as handle, open(output, "w") as out:
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 12:
                continue

            origin = strand_junction(fields, genome, introns)
            junctions[origin] += 1
            reads[origin] += int(fields[4])
            out.write("\\t".join(fields) + "\\n")

    write_summary(sample, junctions, reads, summary)


if __name__ == "__main__":
    prefix = "${task.ext.prefix}" if "${task.ext.prefix}" not in ["null", ""] else "${meta.id}"

    os.makedirs("stranded", exist_ok=True)
    main("${meta.id}", "${junc}", "${fasta}", "${fai}", "${gtf}", f"stranded/{prefix}.junc", f"{prefix}.strand_summary.tsv")

    versions = {"${task.process}": {"python": platform.python_version()}}
    with open("versions.yml", "w") as handle:
        yaml.dump(versions, handle)
