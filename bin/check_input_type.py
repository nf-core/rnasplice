#!/usr/bin/env python3

import argparse
import pandas as pd
import os


# Define command-line arguments
parser = argparse.ArgumentParser(description='Validate a samplesheet file and output files with different column sets')
parser.add_argument('samplesheet_file', help='Path to the samplesheet file')
args = parser.parse_args()

# Define the column sets
fastq_columns = ['sample', 'fastq_1', 'fastq_2', 'strandedness', 'condition']
bam_columns = ['sample', 'bam', 'condition']
transcriptome_columns = ['sample', 'bam', 'transcriptome', 'condition']
salmon_columns = ['sample', 'salmon', 'condition']

# Load the samplesheet file
samplesheet = pd.read_csv(args.samplesheet_file)
base_samplesheet = os.path.splitext(args.samplesheet_file)[0]

# Check that required columns exist and are not empty
required_columns = ['sample', 'condition']
for col in required_columns:
    if col not in samplesheet.columns:
        print(f'Error: required column "{col}" not found in the samplesheet file')
        exit(1)
    if samplesheet[col].isna().any():
        print(f'Error: required column "{col}" has missing values in the samplesheet file')
        exit(1)

# Check which column set to output based on the available columns
if all(col in samplesheet.columns and not samplesheet[col].isna().all() for col in ['fastq_1', 'fastq_2']):
    # Output the fastq column set
    samplesheet[fastq_columns].to_csv(base_samplesheet + '.fastq.csv', index=False)
elif 'bam' in samplesheet.columns and not samplesheet['bam'].isna().all():
    if 'transcriptome' in samplesheet.columns and not samplesheet['transcriptome'].isna().all():
        # Output the transcriptome column set
        samplesheet[transcriptome_columns].to_csv(base_samplesheet + '.transcriptome.csv', index=False)
    else:
        # Output the bam column set
        samplesheet[bam_columns].to_csv(base_samplesheet + '.bam.csv', index=False)
elif 'salmon' in samplesheet.columns and not samplesheet['salmon'].isna().all():
    # Output the salmon column set
    samplesheet[salmon_columns].to_csv(base_samplesheet + '.salmon.csv', index=False)
else:
    print('Error: no valid column set found in the samplesheet file')
    exit(1)

print('Validation successful')
