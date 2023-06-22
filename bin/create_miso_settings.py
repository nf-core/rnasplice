#!/usr/bin/env python
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import argparse
from itertools import repeat
from parsimonious.grammar import Grammar

parser = argparse.ArgumentParser()
parser.add_argument("--bam_prefix", type=str)
parser.add_argument("--miso_prefix", type=str)
parser.add_argument("--bams", nargs="+")
parser.add_argument("--name", type=str, nargs="+")
parser.add_argument("--width", type=int)
parser.add_argument("--height", type=int)
parser.add_argument("--output")

args = parser.parse_args()

bam_prefix = args.bam_prefix
miso_prefix = args.miso_prefix
bams = args.bams
name = args.name
width = args.width
height = args.height
output = args.output

name = "".join(name)
name = name.replace(",", "', '")
name = name.replace("[", "['")
name = name.replace("]", "']")

n = len(bams)
col = "#CC0011"
n_col = list(repeat(col, n))


settings = f"""[data]
bam_prefix = {bam_prefix}/
miso_prefix = {miso_prefix}/
bam_files = {bams}
miso_files = {name}

[plotting]
fig_width = {width}
fig_height = {height}
intron_scale = 30
exon_scale = 4
logged = False
font_size = 6
ymax = 500
show_posteriors = True
bar_posteriors = False
number_junctions = True
resolution = .5
posterior_bins = 40
gene_posterior_ratio = 5
colors = {n_col}
bar_color = "b"
bf_thresholds = [0, 1, 2, 5, 10, 20]
"""

with open(output, "w") as fout:
    fout.writelines(settings)
