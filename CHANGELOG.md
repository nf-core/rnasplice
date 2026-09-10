# nf-core/rnasplice: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

### Added

- Add first steps of leafcutter splicing quantification
- #200 - Add `--ignore_tx_version` to ignore transcript versions in tximport, for tx2gene files whose transcript IDs carry no version (e.g. from GENCODE) (requested by @Oliverfeudj, done by @piplus2)
- #204 - Improve documentation for `--rmats_paired_stats` (requested by @mlbonatelli, done by @piplus2)
- #212 - Add test config for unpaired `rMATS` (by @piplus2)
- #257 - Add a `test_dexseq_gff` profile and a pipeline level nf-test for a user supplied DEXSeq annotation (`--gff_dexseq`) (by @piplus2)
- #258 - Add pipeline level nf-tests for `--source genome_bam`, `transcriptome_bam` and `salmon_results`, which had no test coverage (by @piplus2)
- #261 - Add a pipeline nf-test for samples split over several runs (by @piplus2)
- #263 - Accept optional `strandedness` and `single_end` columns when starting from `--source genome_bam` or `transcriptome_bam`, defaulting to `unstranded` and paired end (requested by @albamasmalavila, done by @piplus2)

### Changed

- #195 - Sync with the nf-core template 4.0.2, migrate samplesheet and contrastsheet validation from Python scripts to nf-schema, move local modules to `TOOL/SUBTOOL/main.nf`, track software versions with `topic: versions` and replace the deprecated `CUSTOM_GETCHROMSIZES` with `SAMTOOLS_FAIDX` (by @piplus2)
- #195 - `PIPELINE_INITIALISATION` handles all four input sources (fastq, genome_bam, transcriptome_bam, salmon_results) (by @piplus2)
- #207 - Update `StageR` to 1.32.0, `HTSeq` to 2.1.2 and `DEXSeq` to 1.56.0 (by @piplus2)
- #210 - Bump `nf-schema` to 2.7.2 and apply static typing to `params` in `main.nf`. The minimum Nextflow version is now 26.04.0 (by @piplus2)
- #212 - Default `--rmats_paired_stats` to `false`. The paired test requires a specific paired design in the samplesheet, so it is now opt in, see the documentation (by @piplus2)
- #215 - Refactor the `dexseq` modules to the nf-core module template (by @piplus2)
- #219 - Rename `MISO_INDEX` and `MISO_RUN` to `MISOPY_INDEX` and `MISOPY_RUN` and refactor them to the nf-core module template (by @piplus2)
- #220 - Remove `parse_miso_index.py`, which misopy does not need (by @piplus2)
- #221 - Refactor `SUBREAD_FLATTENGTF` to the nf-core module template (by @piplus2)
- #225 - `DEXSEQ_DTU` runs DEXSeq in parallel through `BPPARAM`, defaulting to `SerialParam()` (by @piplus2)
- #229 - Use the nf-core `SUPPA` modules and move the `SUPPA` helper modules to `modules/local` (by @piplus2)
- #238, #244, #253 - Sync with the nf-core templates 4.0.3, 4.1.0 and 4.1.0-2 (by @piplus2)
- #242 - Refactor the `misopy` modules to the nf-core module template (by @piplus2)
- #243 - Remove the unused `contrast_check` workflow (by @piplus2)
- #245 - Refactor `TXIMPORT` to the nf-core module template, update `bioconductor-tximeta` 1.8.0 -> 1.28.2 (`tximport` 1.38.2, R 4.0 -> 4.5) and default `--ignore_tx_version` to `FALSE` (by @piplus2)
- #246 - Replace the local `LEAFCUTTER_CLUSTER` module and `bin/leafcutter_cluster_regtools.py` with the nf-core `leafcutter/clusterregtools` module (by @piplus2)
- #247 - Rename `GTF_2_GFF3` to `GFFREAD_GTF2GFF3` and refactor it to the nf-core module template (by @piplus2)
- #248 - Move `GFFREAD_TX2GENE` to `modules/local/gffread/tx2gene` and update `gffread` 0.12.1 -> 0.12.7 (by @piplus2)
- #249 - Refactor `STAGER` to the nf-core module template and expose the `stageR` options through `ext.args` (`--alpha`, `--method`, `--allow_na`, `--only_significant_genes`, `--order`), which were previously ignored. Defaults are unchanged (by @piplus2)
- #250 - Move the `VISUALISE_MISO` subworkflow to the nf-core subworkflow template (by @piplus2)
- #251 - Rename `DRIMSEQ_FILTER` to `DRIMSEQ_DMFILTER` and refactor it to the nf-core module template (by @piplus2)
- #254 - Rename `GTF_GENE_FILTER` to `GTFGENEFILTER` and refactor it to the nf-core module template. `bin/filter_gtf_for_genes_in_genome.py` is now a module template (by @piplus2)
- #255 - Refactor `EDGER_EXON` to the nf-core module template. `bin/run_edger_exon.R` is now a module template, `--n_edger_plot` moves to `ext.args` and the inputs and outputs carry a `meta` map (by @piplus2)
- #256 - Move the `ALIGN_STAR` subworkflow to the nf-core subworkflow template, with nf-tests for both the `STAR_ALIGN` and the `STAR_ALIGN_IGENOMES` paths (by @piplus2)
- #257 - Move the `DEXSEQ_DEU` subworkflow to the nf-core subworkflow template (by @piplus2)
- #257 - `DEXSEQ_DEU` emits the DEXSeq exon count tables sorted by file name, and takes an explicit `prepare_annotation` input instead of reading `params.gff_dexseq` (by @piplus2)
- #261 - The pipeline parses the samplesheet and the contrastsheet once each, in `PIPELINE_INITIALISATION` (by @piplus2)
- #261 - Remove the unused `INPUT_CHECK` subworkflow (by @piplus2)
- #264 - Move the `DRIMSEQ_DEXSEQ_DTU` subworkflow to the nf-core subworkflow template (by @piplus2)

### Fixed

- #195 - Fix Nextflow 26+ compatibility: no deprecated `if/else` blocks in `nextflow.config`, `ext.when` closures in `conf/modules.config` and no `switch/case` statements (by @piplus2)
- #195 - Fix channel reuse in the `SUPPA` subworkflow, `CREATE_BAMLIST` null paths for single condition rMATS runs, `SUBREAD_FLATTENGTF` multi-line versions breaking the YAML and the `STAR_ALIGN` argument renamed to `--quantTranscriptomeSAMoutput` (by @piplus2)
- #198 - Do not pass sorted BAM files to Salmon (reported by @albamasmalavila, fix by @piplus2)
- #199 - Fix the ignored `--clip_r1` and `--clip_r2` arguments of `TRIMGALORE` (reported by @misiti-1864309, fix by @piplus2)
- #201 - Pass the contrasts channel, instead of the samplesheet channel, to `RMATS` and `SUPPA` (by @piplus2)
- #208 - Follow the nf-core standard for conditional output in `publishDir.saveAs` (by @piplus2)
- #140 - Add validation of sample names for R compatibility (by @fhausmann)
- #212 - Fix the sample order of the `RMATS` paired model, which was not taken from the samplesheet (by @piplus2)
- #214 - Fix the missing `BPPARAM` argument of `DEXSEQ_DTU` and remove its unused `ntop` argument (by @piplus2)
- #218 - Fix a regex error in `parse_miso_index.py` (by @piplus2)
- #219 - Fix `MISOPY` handling of paired-end reads and the BAM files passed to it (by @piplus2)
- #224 - Fix the sample IDs in the `DEXSEQ_DTU` script (reported by @albamasmalavila, fix by @piplus2)
- #232 - Fix the wrong branch taken when the input is BAM (reported by @albamasmalavila, fix by @piplus2)
- #238 - Add a missing backslash in the `mergeevents` command (by @piplus2)
- #241 - Add a missing backslash in the `clustergroups` command (by @piplus2)
- #245 - Fix the ignored `--tximport_ignore_tx_version` argument of `TXIMPORT` (by @piplus2)
- #254 - Fix `GTFGENEFILTER` logging the characters of the last GTF sequence name instead of the sequence names found, and raising a `NameError` on an empty GTF. It now fails explicitly when no GTF feature matches the genome FASTA (by @piplus2)
- #255 - Fix `EDGER_EXON` writing an unreadable zero-page PDF when a contrast has no gene to plot. It now fails explicitly on missing samplesheet or contrastsheet columns, on a sample assigned to more than one condition, on a contrast condition absent from the samplesheet and on a missing featureCounts table (by @piplus2)
- #256 - Fix the software version commands of `STAR_ALIGN_IGENOMES` and `STAR_GENOMEGENERATE_IGENOMES`, which failed with a shell syntax error and aborted both processes (by @piplus2)
- #256 - Use `--quantTranscriptomeBan Singleend` on the iGenomes path, as the pinned STAR 2.6.1d does not accept `--quantTranscriptomeSAMoutput` (by @piplus2)
- #257 - Fix `--gff_dexseq` aborting with `Not a valid argument for 'combine' operator [String]`. `RNASPLICE` overwrote the `ch_dexseq_gff` channel it receives from `PREPARE_GENOME` with the raw parameter string (by @piplus2)
- #258 - Fix `--source genome_bam`, `transcriptome_bam` and `salmon_results` aborting at startup, as the samplesheet was always validated against the fastq schema. The schema is now selected from `--source` by `samplesheetSchema()` (by @piplus2)
- #258 - `PIPELINE_INITIALISATION` parses the samplesheet and passes it through `NFCORE_RNASPLICE` into `RNASPLICE`, as in the nf-core template. `RMATS` receives that samplesheet channel instead of a second `file(params.input)`, and the rMATS single condition check goes through the same source schema instead of reading the CSV with `new File`/`new URL`, which only handled local and `http` paths (by @piplus2)
- #258 - Fix `BAM_SORT_STATS_SAMTOOLS` being called with a bare fasta path in the `genome_bam` and `transcriptome_bam` paths, which aborted the run (by @piplus2)
- #258 - Fix `No such variable: ch_txi_suppa_tpm` with `--source salmon_results` and `--suppa` (by @piplus2)
- #258 - Fix `ISOFORMSWITCHANALYZER` receiving the `.tar.gz` archives instead of the extracted Salmon directories with `--source salmon_results`. `TX2GENE_TXIMPORT` extracts them and now emits them (by @piplus2)
- #261 - Fix the paired rMATS sample order when a sample spans several samplesheet rows (by @piplus2)
- #261 - Fix single condition rMATS runs, which aborted before producing any output (by @piplus2)
- #261 - Fix samples split over several samplesheet rows, whose fastq files were passed on without being concatenated (by @piplus2)
- #263 - **Breaking change**: `DEXSEQ_COUNT` no longer counts BAM input as forward stranded, it follows the samplesheet `strandedness`, so DEXSeq results change for BAM samplesheets that do not set the column (reported by @albamasmalavila, fix by @piplus2)
- #264 - Fix the `DEXSEQ_DTU` stub, which wrote file names the process outputs did not match, so a stub run of the DTU path failed (by @piplus2)

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
