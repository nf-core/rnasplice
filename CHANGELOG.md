# nf-core/rnasplice: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

### Added

- Add first steps of leafcutter splicing quantification
- #200 -Ignore transcript version in tximport by default to improve compatibility with tx2gene files from GTFs (e.g., from GENCODE) that may not include version numbers in transcript IDs. This can be overridden with `--ignore_tx_version false` if desired (requested by @Oliverfeudj, done by @piplus2)
- #204 - Improve documentation for `--rmats_paired_stats` (requested by @mlbonatelli, done by @piplus2)
- #212 - Add test config for unpaired `rMATS` (by @piplus2)

### Changed

- #195 -Synced pipeline with nf-core template 4.0.2 (by @piplus2)
- #195 - Migrated samplesheet and contrastsheet validation from Python scripts to nf-schema (by @piplus2)
- #195 - Replaced deprecated `CUSTOM_GETCHROMSIZES` module with `SAMTOOLS_FAIDX` (by @piplus2)
- #195 - Migrated all modules to use `topic: versions` pattern for software version tracking (by @piplus2)
- #195 - Migrated local modules to `TOOL/SUBTOOL/main.nf` directory structure (by @piplus2)
- #195 - Updated `PIPELINE_INITIALISATION` to handle all 4 input source types (fastq, genome_bam, transcriptome_bam, salmon_results) (by @piplus2)
- #207 - Updated `StageR` to version 1.32.0 for improved performance and bug fixes (by @piplus2)
- #207 - Updated `HTSeq` to version 2.1.2 for improved performance and bug fixes (by @piplus2)
- #207 - Updated `DEXSeq` to version 1.56.0 for improved performance and bug fixes (by @piplus2)
- #210 - Bumped `nf-schema` to version 2.7.2 (@piplus2)
- #210 - Applied static typing to `params` in `main.nf` (by @piplus2)
- Minimum Nextflow version updated to 26.04.0 for static typing support.
- #212 - `RMATS` paired stats default changed to `false`, as this is not appropriate for all datasets. The paired differential splicing test requires a specific paired design in the samplesheet, and users should explicitly enable this option if their dataset meets the requirements. See documentation for details. (by @piplus2)
- #215 - Refactored `dexseq` modules to match nf-core schema (by @piplus2)
- Refactored `MISO_INDEX` module into `MISOPY_INDEX` to match nf-core module template.
- #219 Refactored `MISO_RUN` module into `MISOPY_RUN` to match nf-core module template (by @piplus2)
- #220 - Removed `parse_miso_index.py` script, as it is not necessary for misopy to work (by @piplus2)
- #221 - Refactored `SUBREAD_FLATTENGTF` to match nf-core module template (by @piplus2)
- #225 -Now `DEXSEQ_DTU` module uses `BPPARAM` argument to enable parallel processing, if no `BPPARAM` is provided, it defaults to `SerialParam()` (by @piplus2)
- #229 - Use nf-core modules for `SUPPA` and moved `SUPPA` helper modules to `modules/local` using templates (by @piplus2)
- #238 - Synced pipeline with nf-core template 4.0.3 (by @piplus2)
- #242 - Refactored `misopy` modules to match nf-core module template (by @piplus2)
- #243 - Removed unused workflows `contrast_check` (by @piplus2)
- #244 - Updated to nf-core template 4.1.0 (by @piplus2)

### Fixed

- #195 - Fixed Nextflow 26+ compatibility: removed deprecated `if/else` blocks from `nextflow.config` (by @piplus2)
- #195 - Fixed Nextflow 26+ compatibility: migrated `conf/modules.config` from `if` blocks to `ext.when` closures (@piplus2)
- #195 - Fixed Nextflow 26+ compatibility: replaced `switch/case` statements with `if/else if` chains (by @piplus2)
- #195 - Fixed channel reuse bugs in SUPPA subworkflow (by @piplus2)
- #195 - Fixed `CREATE_BAMLIST` null path handling for single-condition rMATS runs (by @piplus2)
- #195 - Fixed `SUBREAD_FLATTENGTF` multi-line version output breaking YAML parsing (by @piplus2)
- #195 - Fixed `STAR_ALIGN` deprecated `--quantTranscriptomeBan` parameter renamed to `--quantTranscriptomeSAMoutput` (by @piplus2)
- #198 - Fixed a bug that passed sorted BAM files to Salmon (reported by @albamasmalavila, fix by @piplus2)
- #199 - Fixed ignored arguments `--clip_r1` and `--clip_r2` in `TRIMGALORE` module for NextSeq trimming and read clipping (reported by @misiti-1864309, fix by @piplus2)
- #201 - Pass the correct contrasts channel instead of the samplesheet channel `RMATS` and `SUPPA` (by @piplus2)
- #208 - Modules `publishDir.saveAS` follow the nf-core standard for conditional output by (by @piplus2)
- #140 - Add validation for sample names to be compliant with R (by @fhausmann)
- #212 - Fixed bug in `RMATS` where the order of the samples in the paired model was not correctly determined from the samplesheet (by @piplus2)
- #214 - Fixed missing `BPPARAM` argument in `DEXSEQ_DTU` module to run DEXSeq with parallel processing (by @piplus2)
- #214 - Removed unused `ntop` argument from `DEXSEQ_DTU` module and script (by @piplus2)
- #218 - Fixed regex error in `parse_miso_index.py` (by @piplus2)
- #219 - Fixed `MISOPY` module to correctly handle paired-end reads and pass the correct BAM files to MISOPY (by @piplus2)
- #224 - Fixed `DEXSEQ_DTU` script to correctly handle the sample IDs (reported by @albamasmalavila, fix by @piplus2)
- #232 - Fixed wrong branch when input is BAM in `dev` (reported by @albamasmalavila, fix by @piplus2)
- #238 - Added missing backslash in `mergeevents` command in `dev` branch (by @piplus2)
- #241 - Added missing backslash in `clustergroups` command in `dev` branch (by @piplus2)

## v1.0.5 - 2024-11-03

- Added IsoformSwitchAnalyzeR to pipeline.
- Updated to migrate from nf-validation to nf-schema.
- Updated to remove max_memory, max_cpus, max_time and replace with process resourceLimits.
- Updated for nf-core template version 3.0.2.
- Updated to remove lib folder and add utils subworkflows.

## Dev v1.0.5dev - TBD - TBD

- Add first steps of leafcutter splicing quantification

## v1.0.4 - 2024-04-21

- Fixed incorrect assignment of cluster groups (Issue #131).

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
