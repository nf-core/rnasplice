#!/usr/bin/env python3
"""
Give a strand to the junctions regtools could not strand.

An unstranded library has no read orientation to take the junction strand from, and
regtools writes `?` in the strand column unless the aligner left an `XS` tag behind.
The strand of a junction is however written in the genome itself: an intron whose
first and last two bases spell one of the canonical splice motifs can only be read in
one direction. Junctions with a non-canonical motif are looked up in the annotation,
and the ones that are in neither are left as they are, for the clustering to discard.
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


def main(junc, fasta, fai, gtf, output):
    genome = IndexedFasta(fasta, fai)
    introns = annotated_introns(gtf)
    counts = {"total": 0, "already stranded": 0, "from the motif": 0, "from the annotation": 0, "left strandless": 0}

    with open(junc) as handle, open(output, "w") as out:
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 12:
                continue
            counts["total"] += 1

            if fields[5] != "?":
                counts["already stranded"] += 1
                out.write(line)
                continue

            bounds = intron_bounds(fields)
            if bounds is None:
                counts["left strandless"] += 1
                out.write(line)
                continue

            start, end = bounds
            motif = (genome.fetch(fields[0], start, start + 1), genome.fetch(fields[0], end - 1, end))
            strand = CANONICAL_MOTIFS.get(motif)
            if strand:
                counts["from the motif"] += 1
            else:
                strand = introns.get((fields[0], start, end))
                if strand:
                    counts["from the annotation"] += 1
                else:
                    counts["left strandless"] += 1
                    out.write(line)
                    continue

            fields[5] = strand
            out.write("\\t".join(fields) + "\\n")

    for name, count in counts.items():
        print(f"{name}: {count}")


if __name__ == "__main__":
    prefix = "${task.ext.prefix}" if "${task.ext.prefix}" not in ["null", ""] else "${meta.id}"

    os.makedirs("stranded", exist_ok=True)
    main("${junc}", "${fasta}", "${fai}", "${gtf}", f"stranded/{prefix}.junc")

    versions = {"${task.process}": {"python": platform.python_version()}}
    with open("versions.yml", "w") as handle:
        yaml.dump(versions, handle)
