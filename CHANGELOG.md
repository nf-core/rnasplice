# nf-core/rnasplice: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.3 - 2024-02-23

- Improved TPM file splitting performance (Issue #120).
- Fixed an issue where R scripts altered sample names upon loading (Issue #122).

## v1.0.2 - 2024-01-08

Patch for run_stager.R (#108) and template update v2.11.1 (#109).

## v1.0.1 - 2023-11-15

Patch for run_drimseq_filter.R to cast command line arguments to numeric. See issue #98 on nf-core/rnasplice.

## v1.0.0 - 2023-05-22

First release of nf-core/rnasplice, created with the [nf-core](https://nf-co.re/) template.

### `Added`

Implemented pipeline:

- Merge re-sequenced FastQ files (cat)
- Read QC (FastQC)
- Adapter and quality trimming (TrimGalore)
- Alignment with STAR:
  - STAR -> Salmon
  - STAR -> featureCounts
  - STAR -> HTSeq (DEXSeq count)
- Sort and index alignments (SAMtools)
- Create bigWig coverage files (BEDTools, bedGraphToBigWig)
- Pseudo-alignment and quantification (Salmon; optional)
- Summarize QC (MultiQC)
- Differential Exon Usage (DEU):
  - HTSeq -> DEXSeq
  - featureCounts -> edgeR
  - Quantification with featureCounts or HTSeq
- Differential exon usage with DEXSeq or edgeR
  - Differential Transcript Usage (DTU):
  - Salmon -> DRIMSeq -> DEXSeq
  - Filtering with DRIMSeq
- Differential transcript usage with DEXSeq
- Event-based splicing analysis:
  - STAR -> rMATS
  - Salmon -> SUPPA2

Updated pipeline:

- Visualization of differential results with edgeR, DEXSeq, and MISO
- Contrasts specified using contrastsheet.csv
- Allow users to specify input data type and start point (e.g., fastq, genome_bam, transcript_bam, salmon_results)
- Pipeline schematic updated
