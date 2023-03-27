#!/usr/bin/env python

import argparse
import csv
import os
import re
import sys


def validate_samplesheet(samplesheet_path):
    required_columns = ["sample", "condition", "bam"]
    optional_columns = ["transcriptome"]
    bam_pattern = re.compile(r"\.bam$")

    with open(samplesheet_path, "r") as samplesheet_file:
        reader = csv.DictReader(samplesheet_file)

        # Check that all required columns are present
        missing_columns = set(required_columns) - set(reader.fieldnames)
        if missing_columns:
            sys.exit(f"Error: Missing required columns: {', '.join(missing_columns)}")

        # Check that required columns are not empty and bam column has the right format
        for row in reader:
            for column in required_columns:
                if not row[column]:
                    sys.exit(f"Error: Required column {column} is empty in row {row}")
            if not bam_pattern.search(row["bam"]):
                sys.exit(f"Error: Invalid BAM file format in row {row}")

            # Check transcriptome column if present
            if "transcriptome" in row:
                if not row["transcriptome"]:
                    sys.exit(f"Error: Transcriptome column is empty in row {row}")
                if not bam_pattern.search(row["transcriptome"]):
                    sys.exit(f"Error: Invalid transcriptome file format in row {row}")

    # Write valid samplesheet file
    valid_samplesheet_path = os.path.splitext(samplesheet_path)[0] + ".valid.csv"
    os.rename(samplesheet_path, valid_samplesheet_path)
    print(f"Valid samplesheet file written to {valid_samplesheet_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate a samplesheet file")
    parser.add_argument("samplesheet", help="Path to the samplesheet file")
    args = parser.parse_args()

    validate_samplesheet(args.samplesheet)
