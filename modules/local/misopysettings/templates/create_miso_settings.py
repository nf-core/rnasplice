#!/usr/bin/env python
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

import platform
from itertools import repeat

import yaml


def write_settings(bam_prefix, miso_prefix, bams, miso_data, width, height, output):
    n = len(bams)
    col = "#CC0011"
    n_col = list(repeat(col, n))

    # NB: Assume bam_files and miso_files have same prefix. Sorting independently
    # should give same sort order for both. Requirement for MISO.

    settings = f"""[data]
bam_prefix = {bam_prefix}/
miso_prefix = {miso_prefix}/
bam_files = {sorted(bams)}
miso_files = {sorted(miso_data)}

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


def main():
    bam_prefix = "."
    miso_prefix = "."
    ext_args_str = "${task.ext.args}"
    if ext_args_str:
        # transform string into dictionary (arg, value)
        ext_args = {}
        # split by space, then loop through each argument and its value
        split_args = ext_args_str.split()
        for i in range(0, len(split_args), 2):
            arg = split_args[i]
            value = split_args[i + 1] if i + 1 < len(split_args) else None
            ext_args[arg] = value

        if "--bam_prefix" in ext_args:
            bam_prefix = ext_args["--bam_prefix"]

        if "--miso_prefix" in ext_args:
            miso_prefix = ext_args["--miso_prefix"]

    list_bams = str("${bams}").split(" ")
    list_miso_data = str("${miso_data}").split(" ")

    write_settings(
        bam_prefix,
        miso_prefix,
        list_bams,
        list_miso_data,
        int("${fig_width}"),
        int("${fig_height}"),
        "miso_settings.txt",
    )

    # Write version information to versions.yml
    versions = {"${task.process}": {"python": platform.python_version()}}

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)


if __name__ == "__main__":
    main()
