#!/usr/bin/env python3
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import argparse
import itertools

def main(args):

    handle = open(args.file, 'r')

    header = handle.readline()

    samples = header.split('\t')

    conditions = [sample.rsplit('_', 1)[0] for sample in samples]

    handle.close()

    # # #

    last_index = 0
    out = []
    for v, g in itertools.groupby(enumerate(conditions), lambda k: k[1]):
        l = [*g]
        out.append([last_index + 1, l[-1][0] + 1])
        last_index += len(l)

    assert len(out) == 2, "Column numbers have to be continuous, with no overlapping or missing columns between them."

    # # #

    groups = [f'{start}-{end}' for start, end in out]

    groups = ','.join(groups)

    print(groups, end = "")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()
    main(args)
