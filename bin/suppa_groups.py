#!/usr/bin/env python3
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import argparse
import itertools


def main(args):
    # Open file
    handle = open(args.file, "r")

    # Read header only
    header = handle.readline()

    # Split header by delimiter (e.g., transcript_GBR_1)
    samples = header.split("\t")

    # Trim replicate number (e.g., transcript_GBR_1 -> transcript_GBR)
    conditions = [sample.rsplit("_", 1)[0] for sample in samples]

    # Close file
    handle.close()

    # Create list of consecutive condition indices (e.g., [[1,2], [3,4]])
    last_index = 0
    out = []
    for v, g in itertools.groupby(enumerate(conditions), lambda k: k[1]):
        l = [*g]
        out.append([last_index + 1, l[-1][0] + 1])
        last_index += len(l)

    # Assert that condition indices are consecutive
    assert len(out) == 2, "Column numbers have to be continuous, with no overlapping or missing columns between them."

    # Format ranges for printing
    groups = ",".join([f"{start}-{end}" for start, end in out])

    # Print to stdout without newline
    print(groups, end="")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()
    main(args)
