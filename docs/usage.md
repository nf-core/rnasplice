# nf-core/rnasplice: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnasplice/usage](https://nf-co.re/rnasplice/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The samplesheet has to be a comma-seperated file with at least 3 columns, and a header row as shown in the example below.

```csv
sample,fastq_1,fastq_2,strandedness,condition
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward,control
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward,control
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward,control
```

The samplesheet needs the `sample` and `condition` columns, and at least one of the input columns `fastq_1`, `genome_bam`, `transcriptome_bam` or `salmon_results`. See [Input types](#input-types) for what each input column is used for.

| Column              | Description                                                                                                                                                                                                                     |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`            | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. It must be a valid R name: letters, digits, `.` and `_`, starting with a letter or `.`.                           |
| `condition`         | The name of the condition a sample belongs to (e.g. 'control', or 'treatment') - these labels will be used for downstream analysis.                                                                                             |
| `fastq_1`           | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                      |
| `fastq_2`           | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                      |
| `strandedness`      | Sample strand-specificity. Must be one of `unstranded`, `forward`, `reverse` or `auto` to infer it, see [Strandedness](#strandedness). Required for samples with FASTQ files. Optional otherwise, defaulting to `unstranded`.   |
| `single_end`        | Read type of a sample without FASTQ files, `true` for single end and `false` for paired end. Optional, defaults to `false`. When `fastq_1` is given, the read type is taken from the presence of `fastq_2` and this is ignored. |
| `genome_bam`        | Full path to a genome BAM file from a splice aware aligner, such as the `*.sorted.bam` files of nf-core/rnaseq. The pipeline sorts and indexes it. File has to be in ".bam" format.                                             |
| `transcriptome_bam` | Full path to a transcriptome BAM file, such as the `*.Aligned.toTranscriptome.out.bam` files of nf-core/rnaseq. File has to be in ".bam" format.                                                                                |
| `salmon_results`    | Full path to a `.tar.gz` archive of the result folder of a Salmon quantification.                                                                                                                                               |

The `seq_platform`, `seq_center` and `percent_mapped` columns that nf-core/rnaseq writes are accepted and ignored, so its samplesheets can be used as they are, see [Starting from an nf-core/rnaseq run](#starting-from-an-nf-corernaseq-run).

In the case of paired differential splicing analysis (e.g., with `rMATS`), the treatment and control samples order must match:

```csv
sample,fastq_1,fastq_2,strandedness,condition
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward,treatment
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward,treatment
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward,treatment
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz,reverse,control
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz,reverse,control
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz,reverse,control
```

### Strandedness

The `strandedness` of a sample can be set to `auto` to have the pipeline infer it. With FASTQ input this is done as in [nf-core/rnaseq](https://nf-co.re/rnaseq): the reads of each `auto` sample are sub-sampled to 1 million with [fq](https://github.com/stjude-rust-labs/fq) and pseudo-aligned with [Salmon](https://salmon.readthedocs.io/en/latest/library_type.html) in its library type detection mode, before alignment. The fraction of reads assigned to each strand then decides the strandedness:

- `forward` or `reverse` when at least `--stranded_threshold` (default `0.8`) of the stranded reads are on that strand.
- `unstranded` when the forward and reverse fractions differ by at most `--unstranded_threshold` (default `0.1`).
- Otherwise the strandedness is undetermined. The sample is analysed as `unstranded`, with a warning in the log.

The inferred value is then used for every step that needs one: Salmon, featureCounts, DEXSeq, rMATS, LeafCutter and ASpli. A `forward`, `reverse` or `unstranded` value in the samplesheet is used as given, so samples set to `auto` and samples with a given strandedness can be mixed in the same samplesheet. Every run of a sample must have the same `strandedness`.

BAM input has no reads to sub-sample, so `auto` is handled differently:

- With genome BAM input, [RSeQC infer_experiment.py](https://rseqc.sourceforge.net/#infer-experiment-py) is run on the BAM file of each `auto` sample, with a BED12 gene model made from the GTF, and its fractions are called with the same thresholds as above. The result is used by DEXSeq, featureCounts (edgeR), rMATS, LeafCutter, ASpli and the bigWig coverage tracks.
- With transcriptome BAM input, the only tool that uses the strandedness is Salmon, which quantifies an `auto` sample with `--libType A` and detects the library type itself. Salmon writes no fragment counts in alignment mode, so the thresholds above do not apply and the strandedness follows Salmon's own call. Setting `--salmon_quant_libtype` overrides this for every sample.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness,condition
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto,control
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,reverse,control
```

The strandedness used can be checked against the alignments with `--strand_check`, which runs [RSeQC infer_experiment.py](https://rseqc.sourceforge.net/#infer-experiment-py) on the genome BAM files with the same thresholds. A sample whose RSeQC result does not match the strandedness it was analysed with is reported with a warning in the log and in the Strandedness Checks table of the MultiQC report. `--strand_check` needs genome alignments, so it works with `--aligner star` or `star_salmon`, and with genome BAM input. A genome BAM sample set to `auto` already has its RSeQC result, so RSeQC is not run on it again, and the table only reports whether a strandedness could be called.

> [!NOTE]
> rMATS needs every sample to have the same strandedness. With `auto` samples this is checked once their strandedness is inferred, and the pipeline stops if the samples differ.

### Input types

Which input columns are filled decides how each analysis gets its input. An input column must be filled for every sample or for none, as the samples of a contrast are compared with each other, and the pipeline stops if it is filled for only some of them.

| Input               | Analyses                                                                                                     |
| ------------------- | ------------------------------------------------------------------------------------------------------------ |
| `fastq_1`           | All, through STAR alignment and Salmon quantification, as configured with `--aligner` and `--pseudo_aligner` |
| `genome_bam`        | DEXSeq exon usage, edgeR, rMATS, LeafCutter, ASpli, MISO sashimi plots and bigWig coverage tracks            |
| `transcriptome_bam` | Salmon quantification in alignment mode, then DEXSeq DTU, SUPPA and IsoformSwitchAnalyzeR                    |
| `salmon_results`    | DEXSeq DTU, SUPPA and IsoformSwitchAnalyzeR                                                                  |

A BAM file or Salmon result given in the samplesheet is used instead of making it from the FASTQ files:

- With `genome_bam`, STAR is not run, and the genome-based analyses use the given BAM files.
- With `transcriptome_bam` or `salmon_results`, the FASTQ files are not pseudo-aligned with Salmon, and the transcript-based analyses use the given files.
- The FASTQ files are only trimmed and quality checked when an analysis still needs them.

`transcriptome_bam` and `salmon_results` cannot be given together, as either one is the input of the transcript-based analyses.

Several rows of the same sample are merged: their FASTQ files are concatenated, and their BAM files and Salmon results must be the same on every row.

The BAM files are checked against the references the pipeline uses: the reference sequences named in the header of each genome BAM file must be in `--fasta`, and those of each transcriptome BAM file in the transcript FASTA, with the same lengths. The pipeline stops with an error naming the sequences that do not match.

For example, genome and transcriptome BAM files together:

```csv
sample,condition,strandedness,genome_bam,transcriptome_bam
CONTROL_REP1,control,reverse,AEG588A1.sorted.bam,AEG588A1.Aligned.toTranscriptome.out.bam
CONTROL_REP2,control,reverse,AEG588A2.sorted.bam,AEG588A2.Aligned.toTranscriptome.out.bam
CONTROL_REP3,control,reverse,AEG588A3.sorted.bam,AEG588A3.Aligned.toTranscriptome.out.bam
TREATMENT_REP1,treatment,reverse,AEG588A4.sorted.bam,AEG588A4.Aligned.toTranscriptome.out.bam
TREATMENT_REP2,treatment,reverse,AEG588A5.sorted.bam,AEG588A5.Aligned.toTranscriptome.out.bam
TREATMENT_REP3,treatment,reverse,AEG588A6.sorted.bam,AEG588A6.Aligned.toTranscriptome.out.bam
```

A BAM file does not record the read type of the library, so without FASTQ files it comes from the `single_end` column, which defaults to paired end. Set it to `true` for single end libraries.

Salmon results:

```csv
sample,condition,salmon_results
CONTROL_REP1,control,AEG588A1.tar.gz
CONTROL_REP2,control,AEG588A2.tar.gz
CONTROL_REP3,control,AEG588A3.tar.gz
TREATMENT_REP1,treatment,AEG588A4.tar.gz
TREATMENT_REP2,treatment,AEG588A5.tar.gz
TREATMENT_REP3,treatment,AEG588A6.tar.gz
```

> [!NOTE]
> The `--source` parameter is deprecated and ignored, as the input type is worked out from the samplesheet columns. It will be removed in a future release.

### Starting from an nf-core/rnaseq run

nf-core/rnaseq writes a samplesheet listing the BAM files of each sample when it is run with `--save_align_intermeds` and the `star_salmon` aligner:

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir $PWD/rnaseq_results \
    --aligner star_salmon \
    --save_align_intermeds \
    --save_reference \
    -profile docker
```

Add a `condition` column to `rnaseq_results/samplesheets/samplesheet_with_bams.csv` and give it to rnasplice with the references nf-core/rnaseq saved to `rnaseq_results/genome/` with `--save_reference`: its genome FASTA, GTF and transcript FASTA. nf-core/rnaseq writes the BAM paths relative to its `--outdir`, so give it an absolute path as above, or run rnasplice from the directory nf-core/rnaseq was run from.

```bash
nextflow run nf-core/rnasplice \
    --input samplesheet_with_bams.csv \
    --contrasts contrastsheet.csv \
    --fasta rnaseq_results/genome/<genome FASTA> \
    --gtf rnaseq_results/genome/<GTF> \
    --transcript_fasta rnaseq_results/genome/genome.transcripts.fa \
    --outdir rnasplice_results \
    -profile docker
```

The genome BAM files then go to the genome-based analyses and the transcriptome BAM files through Salmon to the transcript-based ones, in one run, with no second alignment. The strandedness in the samplesheet is the one nf-core/rnaseq inferred.

The transcriptome BAM files are only valid against the exact transcript FASTA nf-core/rnaseq aligned them to, which `--save_reference` writes to `rnaseq_results/genome/`. Without `--transcript_fasta`, rnasplice builds its own from the genome and the GTF, which can differ, for example when nf-core/rnaseq filtered the GTF. The pipeline then stops at the BAM check above rather than quantify against the wrong transcripts. The saved genome FASTA and GTF are the ones the BAM files were aligned to, which can differ from the ones given to nf-core/rnaseq, for example with `--additional_fasta`, which adds sequences such as ERCC spike-ins to both, and to the BAM headers. The BAM check stops the run against the original genome FASTA then.

## Contrastsheet input

You will also need to create a contrastsheet with information about the contrasts you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--contrasts '[path to contrastsheet file]'
```

The contrastsheet has to be a comma-separated file with a header row, using the same column names as [nf-core/differentialabundance](https://nf-co.re/differentialabundance):

```csv
id,variable,reference,target
TREATMENT_CONTROL,condition,CONTROL,TREATMENT
```

| Column      | Description                                                                                                          |
| ----------- | -------------------------------------------------------------------------------------------------------------------- |
| `id`        | An arbitrary identifier, unique within the contrastsheet, used to name contrast-wise output files.                   |
| `variable`  | The samplesheet column to compare. Only `condition` is supported.                                                    |
| `reference` | The reference (base) level of the comparison, a `condition` of the samplesheet.                                      |
| `target`    | The target level of the comparison, a `condition` of the samplesheet. Log fold changes are target over reference.    |
| `blocking`  | Optional. Samplesheet columns to adjust for, separated by semicolons. See [Blocking variables](#blocking-variables). |

The contrastsheet can have more columns, which are ignored.

### Blocking variables

The `blocking` column names samplesheet columns to adjust for in a contrast, such as batch or donor, separated by semicolons without spaces. Each contrast can have its own blocking variables, or none:

```csv
id,variable,reference,target,blocking
TREATMENT_CONTROL,condition,CONTROL,TREATMENT,batch;donor
```

Only some methods can adjust for blocking variables. The others ignore them, and the pipeline prints a warning when the contrastsheet sets any, so take this into account when you compare results between methods:

| Method                                            | Adjusts for blocking variables |
| ------------------------------------------------- | ------------------------------ |
| DEXSeq DEU (`--dexseq_exon`)                      | Yes                            |
| DEXSeq DTU (`--dexseq_dtu`)                       | Yes                            |
| edgeR DEU (`--edger_exon`)                        | Yes                            |
| rMATS (`--rmats`)                                 | No                             |
| SUPPA (`--suppa`)                                 | No                             |
| IsoformSwitchAnalyzeR (`--isoformswitchanalyzer`) | No                             |
| ASpli (`--aspli`)                                 | No                             |
| LeafCutter (`--leafcutter`)                       | No                             |

DEXSeq adds each blocking variable to its models as an interaction with exon, for example `~ sample + exon + batch:exon + condition:exon` against `~ sample + exon + batch:exon`. The DEXSeq log fold changes are not adjusted for the blocking variables. edgeR adds them to its design matrix, `~ 0 + condition + batch`.

Blocking variables are always treated as categorical, so numeric donor IDs are not read as continuous values. Each blocking column must have a value for every sample, and the same value in every row of a sample with several sequencing runs. They are checked against the samples of each contrast:

- A variable with a single value in the contrast is left out of the model, with a warning.
- Variables confounded with `condition`, or leaving no residual degrees of freedom, stop the pipeline with an error.

### Deprecated column names

The previous column names `contrast`, `treatment` and `control` are still accepted, but are deprecated and the pipeline prints a warning when it reads them. They map to the new columns as follows:

| Deprecated column | New column  |
| ----------------- | ----------- |
| `contrast`        | `id`        |
| `treatment`       | `target`    |
| `control`         | `reference` |

```csv
contrast,treatment,control
TREATMENT_CONTROL,TREATMENT,CONTROL
```

The two sets of column names cannot be mixed in one contrastsheet.

An [example contrastsheet](https://raw.githubusercontent.com/nf-core/test-datasets/rnasplice/chrX_subset/samplesheet/contrastsheet_differentialabundance.csv) is available in the test datasets.

## Alignment options

The pipeline offers the use of [STAR](https://github.com/alexdobin/STAR) (i.e. `--aligner star`) to map raw FastQ reads to a reference genome and to project the alignments onto the transcriptome. Downstream quantification can also be performed following [STAR](https://github.com/alexdobin/STAR) alignment with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) when `--aligner star_salmon` is chosen. Users may wish to perform [STAR](https://github.com/alexdobin/STAR) alignment without [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) quantification when certain downstream tools (e.g. [rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md)) only require BAM files as input. By default intermediate SAM alignment files are not saved for space efficiency reasons. This behaviour can be overriden with the `--save_align_intermeds` parameter. Users should take note that STAR requires a lot of memory (~38GB Human GRch37).

Although [STAR](https://github.com/alexdobin/STAR) is fast, users may wish to choose an even faster option by choosing to psuedo-align and quantify with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) alone by providing the `--pseudo_aligner salmon` parameter. This may be a good option if users do not need to generate BAM files and downstream tools are compatible with this process (e.g. Differential transcript usage using [DEXSeq](https://f1000research.com/articles/7-952), or Event-based approach [SUPPA](https://github.com/comprna/SUPPA)). It should be noted, however, that by defining the parameter `--pseudo_aligner salmon`, users will run [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) in addition to the standard workflow defined by the `--aligner` parameter. This is in keeping with the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) approach, and Salmon can be run in isolation with the addition of the `--skip_alignment` parameter. Furthermore, the pipeline, like [nf-core/rnaseq](https://github.com/nf-core/rnaseq) will generate transcript fasta and Salmon index files from gtf and genome fasta files by default. Users can supply these files, however, if they wish to override this auto-generation using the `--transcript_fasta` and `--salmon_index` parameters. Additional Salmon paramaters can be specified at run-time if users wish. For example, library preparation protocol (library type) by default is inferred from samplesheet information, but can also be specified using the `--salmon_quant_libtype` parameter. You can find additional library type details in the [Salmon documentation](https://salmon.readthedocs.io/en/latest/library_type.html).

### Maximum intron length

By default STAR aligns across introns up to about 590 kb long, `(2^winBinNbits)*winAnchorDistNbins` with its default settings, whatever the genome. In organisms with short introns, such as _Arabidopsis thaliana_, where introns longer than a few kb are rare, this lets STAR report spurious long introns, which then reach the splicing analyses as junctions. Set `--max_intron_length` to a length suited to the organism:

```bash
--max_intron_length 6000
```

It is passed to STAR as `--alignIntronMax` and `--alignMatesGapMax`, and is the intron limit of LeafCutter (`--leafcutter_max_intron_len`) and ASpli (`--aspli_max_intron_len`) unless those are set. With genome BAM input STAR is not run, and it only sets the LeafCutter and ASpli limits. When it is not set, STAR keeps its own limits, LeafCutter uses 500 kb and ASpli 50 kb.

## Quantification options

As discussed above users may wish to align with STAR by providing the `--aligner` parameter, or pseudo-align with Salmon by providing the `--pseudo_aligner salmon` parameter.

There are 3 methods for quantification after STAR alignment:

Already discussed above is quantification using Salmon by providing the `--aligner star_salmon` parameter. This enables access to tools which require Salmon for quantification (e.g. Differential transcript usage using [DEXSeq](https://f1000research.com/articles/7-952), or Event-based approach [SUPPA](https://github.com/comprna/SUPPA))

The pipeline also enables quantification using [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889). This is needed for differential exon usage analysis using edgeR and is activated when the parameter `--edger_exon` is enabled. Please note that as this is aimed at differential exon usage feature type is set as `exon` and cannot be changed. Please take care to use a suitable attribute to categorize the featureCounts attribute type in your GTF using the option `--gtf_group_features` (default: `gene_id`).

The final quantification method following STAR alignment is with [HTSeq](https://htseq.readthedocs.io/en/master/) which is implemented as part of the [DEXSeq](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#3_Counting_reads) package. This is needed for differential exon usage analysis using DEXSeq and is activated when the parameter `--dexseq_exon` is enabled. Using the `--aggregation` parameter the pipeline will combine overlapping genes into a single aggregate gene. This approach can alternatively be skipped and any exons that overlap other exons from different genes will be skipped. Other important options to take note of are the `--alignment_quality` parameter which can be set by the user and defines the minimum alignment quality required for reads to be included (defined in 5th column of a given SAM file) (default: 10). Prior to quantification, DEXSeq provides an annotation preparation script which takes a GTF file as input and returns a GFF file. Users may instead wish to define their own GFF file and skip this annotation preparation skip by supplying it using the `--gff_dexseq` parameter.

## Reference genome files

Like [nf-core/rnaseq](https://github.com/nf-core/rnaseq) the minimum reference genome requirements are a FASTA and GTF file. All other files required to run the pipeline can be generated from these files. If you wish to build new indices and save the output then the `--save_reference` parameter is required. It should be noted, however, that the sequence files and indices of many common species are available for download from [AWS iGenomes](https://nf-co.re/usage/reference_genomes)). Local copies of reference files may also be provided on the command line or via config files (e.g. `--salmon_index '/path/to/salmon_index.tar.gz'`). Where reference files are compressed (e.g. standard files with the `.gz` extension and indices folders with the `tar.gz` extension) these will be automatically uncompressed prior to use.

For nf-core/rnasplice to run it requires a FASTA file and GTF file. These can be specified manually, however users may wish to use the more simple [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) approach:

- Using the `--genome` parameter followed by the genome build (e.g. GRCh37) will mean the FASTA and GTF file will be pulled from AWS-iGenomes, along with available indices required for a given analysis. Furthermore, a local download of AWS genomes from a previous run can also be used by changing the igenomes path with the `--igenomes_base` parameter.

If a GTF is not available a GFF may be used by specifying the `--gff` parameter. This will convert the GFF file into a GTF.

The STAR indices of AWS iGenomes were built with STAR 2.6.x, whose `genomeParameters.txt` STAR 2.7.4a and later refuse to read. A supplied STAR index is passed through `GENOMEPARAMSUPGRADE` first, which rewrites that metadata and leaves the binary index files alone, so the pipeline runs one STAR version throughout. An index built with a recent STAR is left untouched.

Please note if you are using [GENCODE](https://www.gencodegenes.org/) reference genome files please specify the `--gencode` parameter. This is because reference files which come from GENCODE are different to ENSEMBL reference files and this can impact the running of the pipeline. Specifying this parameter can help to mitigate these differences. Furthermore it should be noted that when using GENCODE reference files if you are running Salmon, the `--gencode` flag will also be passed to the index building step (see [this issue](https://github.com/COMBINE-lab/salmon/issues/15)).

## Differential Exon Usage

Differential exon usage (**DEU**) can be performed with two different branches of the pipeline:

### DEXSeq

Following HTSeq quantification you can estimate the differential exon usage using [DEXSeq](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html). This can be chosen by specifying the `--aligner star` or `--aligner star_salmon` followed by the `--dexseq_exon` parameter.

### edgeR

Following featureCounts quantification differential exon usage can also be performed with [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). This can be chosen by specifying the `--aligner star` or `--aligner star_salmon` followed by the `--edger_exon` parameter.

## Differential Transcript Usage

### DEXSeq

Differential transcript Usage (**DTU**) has been implemented from the workflow by [Love et al., 2018](https://f1000research.com/articles/7-952). This analysis can be chosen by specifying the `--aligner star_salmon` or `--pseudo_aligner salmon` followed by the `--dexseq_dtu` parameter. The parameter `--dtu_txi` should either be set as `--dtu_txi scaledTPM` or `--dtu_txi dtuScaledTPM` for successful **DTU** analysis.

Prior to DEXSeq DTU analysis, filtering of genes and features with low expression is completed using [DRIMSeq](https://rdrr.io/bioc/DRIMSeq/man/dmFilter.html) which comes with a number of parameters which should be set by the user. By default these are all set to 0. For example, `--min_samps_gene_expr` defines the minimal number of samples where genes are required to be expressed, and `--min_gene_expr` defines the minimal level of gene expression for genes to be included in the downstream analysis. Similarly, `--min_samps_feature_expr` and `--min_samps_feature_prop` defines the minimal number of samples where features are required to be expressed at a minimal expression or proportion. This minimum level is further defined by additional filtering parameters `--min_feature_expr` and `--min_feature_prop` respectively. Further details of this filtering process can be see within the **DTU** workflow [here](https://f1000research.com/articles/7-952).

### IsoformSwitchAnalyzeR

IsoformSwitchAnalyzeR as described by [Vitting-Seerup et al., 2019](https://doi.org/10.1093/bioinformatics/btz247) performs a genome-wide analysis of splicing patterns and identifies isoform switches. It reads Salmon quantifications, so it runs with FASTQ input quantified with Salmon, with transcriptome BAM input, or with Salmon results input. Two parameters can be set:
`--isoformswitchanalyzer_alpha` sets the significance threshold and `--isoformswitchanalyzer_dIF` sets the minimum absolute difference in isoform usage.

## Event based approaches

Two predominant event-based approaches have been implemented ([rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md) and [SUPPA](https://github.com/comprna/SUPPA)) in this pipeline and can be accessed through `--aligner` or `--pseudo_aligner` options:

### rMATS

[rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md) (Replicate Multivariate Analysis of Transcript Splicing) is used to detec differential alternative splicing from replicated RNA-Seq data. If `--aligner` is set to `--aligner star` or `--aligner star_salmon` then rMATS can be executed with the `--rmats` parameter.

At current, however, there are some restrictions to running an rMATS analysis:

- Samples need to have the same strandedness and read type (single-end/paired-end).

- The `--rmats_paired_stats` can be set to `true` only if there are two conditions. Samples in samplesheet must follow the paired design order (see [Samplesheet input section](#samplesheet-input) above).

rMATS uses the read length to compute the effective length of each event, which the inclusion levels and the test depend on. When `--rmats_read_len` is not set, it is the median length of the reads of all the samples, taken from the samtools stats of the genome BAM files, as rMATS recommends the median for reads of varying length. Reads of other lengths are still counted. The value used is printed in the pipeline log. Set `--rmats_read_len` to use a different value, for example the length before trimming.

### SUPPA2

You can run [SUPPA](https://github.com/comprna/SUPPA) for analyzing splicing events across conditions following pseudo alignment and quantification with Salmon `--pseudo_aligner salmon` or after STAR alignment and Salmon quantification `--aligner star_salmon`, and when the `--suppa` parameter is supplied.

There are two main options for running an analysis with SUPPA - `--suppa_per_local_event` and `--suppa_per_isoform` (the latter is a DTU approach). When `--suppa_per_local_event` is set to `true`, local AS events are calculated and analyzed. When `--suppa_per_isoform` is set to `true`, transcript isoform events are calculated and analyzed.

#### Event Calculation

Events are calculated from user specified annotation files (e.g. GTF files). The parameter `--pool_genes` should be specified when creating ioe/ioi from annotations that are not loci-based. It should be noted that SUPPA advises users to use Ensembl and Gencode annotations to reduce errors at this stage of the analysis.

The `--local_events` parameter requires users to choose the type of events to focus analysis on. They can be any (or all) from the following list (e.g. `--local_events SE SS MX RI FL`):

- SE: Skipping exon (SE) events
- SS: Alternative 5' (A5) and 3' (A3) splice sites (it generates both)
- MX: Mutually Exclusive (MX) exons
- RI: Retained intron (RI)
- FL: Alternative first (AF) and last (AL) exons (it generates both)

#### PSI Calculation

For `local events`, SUPPA reads the `ioe` file generated in the event calculation step and a transcript expression file with the transcript abundances (TPM units) to calculate the relative abundance (PSI) value per sample for each event. It generates a psi file.

For `transcript isoform events`, SUPPA reads the annotation file and a transcript expression file with the transcript abundances (TPM units) to calculate the relative abundance (PSI) value per sample for each event. It generates a psi file.

#### Differential Splicing Analysis

The `PSI` files and `TPM` files are split based on the condition specified in metadata. e.g., condition1.psi, condition2.psi, condition1.tpm, condition2.tpm.

SUPPA then reads the `PSI` for the events and the transcript expression values from multiple samples, grouped by condition, and the `ioe`/`ioi` file, to calculate the events that are differentially spliced between a pair of conditions.

#### Cluster Analysis

Using `dpsi` file and `psivec` file, events are clustered according to `PSI` values across conditions.

### LeafCutter

[LeafCutter](https://davidaknowles.github.io/leafcutter/) quantifies splicing by clustering the introns that spliced reads span. It runs when `--leafcutter` is set, with `--aligner star` or `--aligner star_salmon`, and with genome BAM input.

Junction extraction (regtools `-M`) and clustering keep introns up to `--leafcutter_max_intron_len` bases long. When it is not set, it is [`--max_intron_length`](#maximum-intron-length), or 500,000 if that is not set either. The clustering keeps clusters supported by at least `--leafcutter_min_cluster_reads` reads (default 10), and drops junctions that carry less than `--leafcutter_min_cluster_ratio` of their cluster's reads (default 0.001). These are `leafcutter-cluster`'s `-l`, `-m` and `-p`, whose own defaults are 100,000, 30 and 0.001.

The clusters are then tested for differential splicing with `leafcutter-ds`, once per contrast in the contrastsheet, with the `control` condition as the baseline so that ΔPSI is treatment minus control. Clusters are labelled with the genes whose exons they overlap in the GTF. The `--leafcutter_ds_*` parameters default to `leafcutter-ds`'s own settings.

`leafcutter-ds` stops with an error if either condition of a contrast has fewer samples than `--leafcutter_ds_min_samples_per_intron` or `--leafcutter_ds_min_samples_per_group`. The pipeline leaves such a contrast out of the differential splicing with a warning, and still publishes the clusters. LeafCutter's authors have checked the calibration of its p-values down to four samples per condition, so be careful when lowering these settings.

Every junction LeafCutter clusters needs a strand, and one that has none is discarded. Where the strand comes from depends on the samplesheet `strandedness`:

- `forward` and `reverse`: from the orientation of the reads, which regtools reads off the BAM file itself. Nothing else is needed, so BAM files given in the `genome_bam` column work as they are.

- `unstranded`: from the splice motif. The two bases at each end of every strandless intron are read out of the genome, and the intron takes `+` for GT-AG, GC-AG and AT-AC and `-` for their reverse complements. An intron whose motif is non-canonical is looked up in the GTF and takes the strand of the transcript it belongs to. An intron that is in neither keeps no strand, and the clustering drops it.

> [!IMPORTANT]
> For **unstranded** libraries, the LeafCutter results of this pipeline are **not identical** to the ones you get by following LeafCutter's own instructions.
>
> LeafCutter documents aligning with STAR's `--outSAMstrandField intronMotif`, which makes STAR write the strand it infers from the splice motif into an `XS` tag on each spliced alignment. This pipeline deliberately does not pass that flag, and infers the strand from the genome afterwards instead. The reason is that the flag does more than add a tag: STAR also **suppresses the spliced alignments whose strand it cannot infer**, and those same BAM files are what rMATS, DEXSeq, edgeR and MISO read. Measured on the pipeline's test dataset, enabling it removes about 0.1% of the alignments and moves about 44 read ends per sample from multi-mapping to uniquely mapped, which crosses the MAPQ threshold `DEXSEQ_COUNT` applies and the multi-mapper filter featureCounts applies. An alternative splicing pipeline should not change the input of four of its analyses because a fifth one was switched on.
>
> What this costs LeafCutter, measured on the same dataset: of the junctions the STAR route strands, this route strands 98.7%. The difference is about 37 junctions per sample, every one of them non-canonical, absent from the annotation and supported by a single read. STAR gives those a strand only because they happen to share an alignment with a canonical junction, so the strand is inherited from a neighbour rather than observed. Junctions that are canonical, or annotated, are stranded identically by both routes.
>
> If you need LeafCutter's documented behaviour exactly, pass the flag yourself with `-c`, keeping the rest of the `STAR_ALIGN` arguments as `conf/modules.config` sets them, and be aware that the alignments the other analyses see change with it:
>
> ```groovy
> process {
>     withName: 'STAR_ALIGN' {
>         ext.args = { [
>             '--quantMode TranscriptomeSAM',
>             '--twopassMode Basic',
>             '--outSAMtype BAM Unsorted',
>             '--readFilesCommand gunzip -c',
>             '--runRNGseed 0',
>             '--outFilterMultimapNmax 20',
>             '--alignSJDBoverhangMin 1',
>             '--outSAMattributes NH HI AS NM MD',
>             '--quantTranscriptomeSAMoutput BanSingleEnd',
>             '--outSAMstrandField intronMotif',
>             params.save_unaligned ? '--outReadsUnmapped Fastx' : ''
>         ].join(' ').trim() }
>     }
> }
> ```

A BAM file given in the `genome_bam` column that already carries `XS` tags is stranded from those tags first, and only the junctions left over go through the motif and the annotation.

## Differential Bin and Junction Usage

### ASpli

[ASpli](https://bioconductor.org/packages/ASpli) runs when `--aspli` is set, with `--aligner star` or `--aligner star_salmon`, or with genome BAM input. It is off by default. The pipeline stops at the start if `--aspli` is set without a route that produces genome BAM files.

ASpli splits each gene of the GTF into bins: exonic and intronic pieces bounded by every annotated exon boundary. It counts reads over genes, bins and junctions in one pass over all the samples, then tests each contrast of the contrastsheet. The differences from the other analyses are:

- Compared with DEXSeq and edgeR, which test exon counts only, ASpli also tests intronic bins and the junctions themselves, so intron retention and unannotated splice sites are covered.
- Compared with rMATS, which tests predefined event types between two groups of samples, ASpli tests bins and junctions with one model fitted to all the samples, the same way as edgeR, and reports junction based PSI and PIR (percent intron retention) values alongside.
- `integrateSignals` combines the bin, junction (`localebased`) and anchor junction (`anchorbased`) evidence for each region, so a region supported by more than one signal stands out.

Enable it when you want splicing changes supported by both bin and junction evidence, or when intron retention matters to the analysis. On a large dataset expect it to take longer and use more memory than DEXSeq or edgeR, since ASpli reads every BAM file into memory in full.

Each contrast becomes a contrast vector over the conditions in the order ASpli reads them: `1` for the treatment, `-1` for the control and `0` for every other condition. The model is fitted to all the samples of the samplesheet, as with edgeR, and log fold changes are treatment over control. With more than 24 conditions, ASpli estimates the junction dispersions from the contrasted conditions plus a random sample of the others, with a fixed seed.

Some things to be aware of:

- All the samples must be single end, or all paired end. The pipeline stops otherwise.
- ASpli counts reads over bins and junctions whatever their strand. The samplesheet `strandedness` only decides how the two reads of a pair are given a strand (`strandMode` 1 for `forward`, 2 for `reverse`, 0 for `unstranded` or a mix).
- Junctions longer than `--aspli_max_intron_len` bases are not counted (`maxISize` of `gbCounts`). When it is not set, it is [`--max_intron_length`](#maximum-intron-length), or else 50,000, the value used in the ASpli vignette, as `gbCounts` has no default for it. The value used is written to `aspli/counts/aspli.settings.tsv`.
- ASpli only counts a read over an exon-intron boundary when it overlaps the boundary by at least the minimum read length. `--aspli_read_len` sets it. When it is not set, it is the 5th percentile of the lengths of all the primary mapped reads of each BAM file, taken from the BAM file with the shortest reads, so that most trimmed reads of varying length are still counted. The value used is written to `aspli/counts/aspli.settings.tsv`.
- The options of the ASpli functions are set with `ext.args` (see [custom configuration](#custom-configuration)). `ASPLI_COUNTS` accepts `--min_anchor` and `--junction_threshold`. `ASPLI_DU` accepts `--min_gen_reads`, `--min_bin_reads`, `--min_rds` and `--min_avg_counts` for the filters, `--strong_filter` for `jDUreport`, and `--bin_fc`, `--bin_fdr`, `--bin_inclusion`, `--bjs_inclusion`, `--bjs_fdr`, `--a_inclusion`, `--a_fdr`, `--l_inclusion` and `--l_fdr` for `integrateSignals`. The defaults are the ASpli ones. This includes `--bjs_inclusion 10.3`, as documented for `integrateSignals`, although the other inclusion filters are fractions between 0 and 1. `ASPLI_COUNTS` also accepts `--max_intron_size`, which the pipeline sets from `--aspli_max_intron_len`.
- With the strong filter (the ASpli default), an event is dropped as soon as one of its junctions fails the filter. On a small dataset this can leave no junction to test, and `ASPLI_DU` then stops with a message suggesting `--strong_filter false`:

```groovy
process {
    withName: 'ASPLI_DU' {
        ext.args = '--strong_filter false'
    }
}
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/rnasplice --input samplesheet.csv --contrasts contrastsheet.csv --outdir <OUTDIR> --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/rnasplice -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
contrasts: './contrastsheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rnasplice
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnasplice releases page](https://github.com/nf-core/rnasplice/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASPLICE:RNASPLICE:FASTQ_ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASPLICE:RNASPLICE:FASTQ_ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASPLICE:RNASPLICE:FASTQ_ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASPLICE:RNASPLICE:FASTQ_ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

:::note
We specify the full process name i.e. `NFCORE_RNASPLICE:RNASPLICE:FASTQ_ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.

If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.
:::

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:
   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

:::note
If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.
:::

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
