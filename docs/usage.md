# nf-core/rnasplice: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnasplice/usage](https://nf-co.re/rnasplice/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2,strandedness,condition
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,unstranded,control
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,unstranded,control
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,unstranded,control
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 5 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice and will therefore be automatically concatenated before further downstream analysis.

```console
sample,fastq_1,fastq_2,strandedness,condition
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward,control
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward,control
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward,control
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,reverse,treatment
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,reverse,treatment
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,reverse,treatment
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,,reverse,treatment
```

| Column         | Description                                                                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `strandedness` | Sample strand-specificity. Must be one of `unstranded`, `forward` or `reverse`.                                                                                                        |
| `condition`    | The name of the condition a sample belongs to (e.g. 'control', or 'treatment') - these labels will be used for downstream analysis.                                                    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Alignment options

The pipeline offers the use of [STAR](https://github.com/alexdobin/STAR) (i.e. `--aligner star`) to map raw FastQ reads to a reference genome and to project the alignments onto the transcriptome. Downstream quantification can also be performed following [STAR](https://github.com/alexdobin/STAR) alignment with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) when `--aligner star_salmon` is chosen. Users may which to perform [STAR](https://github.com/alexdobin/STAR) alignment without [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) quantification when certain downstream tools (e.g. [rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md)) only require BAM files as input. By default intermediate SAM alignment files are not saved for space efficiency reasons. This behaviour can be overriden with the `--save_align_intermeds` parameter. Users should take note that STAR requires a lot of memory (~38GB Human GRch37).

Although [STAR](https://github.com/alexdobin/STAR) is fast, users may wish to choose an even faster option by choosing to psuedo-align and quantify with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) alone by providing the `--pseudo_aligner salmon` parameter. This may be a good option if users do not need to generate BAM files and downstream tools are compatible with this process (e.g. Differential transcript usage using [DEXSeq](https://f1000research.com/articles/7-952), or Event-based approach [SUPPA](https://github.com/comprna/SUPPA)). It should be noted, however, that by defining the parameter `--pseudo_aligner salmon`, users will run [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) in addition to the standard workflow defined by the `--aligner` parameter. This is in keeping with the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) approach, and Salmon can be run in isolation with the addition of the `--skip_alignment` parameter. Furthermore, the pipeline, like [nf-core/rnaseq](https://github.com/nf-core/rnaseq) will generate transcript fasta and Salmon index files from gtf and genome fasta files by default. Users can supply these files, however, if they wish to override this auto-generation using using the `--transcript_fasta` and `--salmon_index` parameters. Additional Salmon paramaters can be specified at run-time if users wish. For example, library preparation protocol (library type) by default is inferred from samplesheet information, but can also be specified using the `--salmon_quant_libtype` parameter. You can find additional library type details in the [Salmon documentation](https://salmon.readthedocs.io/en/latest/library_type.html).

## Quantification options

As discussed above users may wish align with STAR by providing the `--aligner` parameter, or pseudo-align with Salmon by providing the `--pseudo_aligner salmon` parameter.

There are 3 methods for quantification after STAR alignment:

Already discussed above is quantification using Salmon by providing the `--aligner star_salmon` parameter. This enables access to tools which required Salmon for quantification (e.g. Differential transcript usage using [DEXSeq](https://f1000research.com/articles/7-952), or Event-based approach [SUPPA](https://github.com/comprna/SUPPA))

The pipeline also enables quantification using [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889). This is associated with differential exon usage analysis using edgeR and is activated when the parameter `--edger_exon` is enabled. Please note that as this is aimed at differential exon usage feature type is set as `exon` and cannot be changed. Please take care to use a suitable attribute to categorize the featureCounts attribute type in your GTF using the option `--gtf_group_features` (default: `gene_id`).

The final quantification method following STAR alignment is with DEXSeq count python script which uses [HTSeq](https://htseq.readthedocs.io/en/master/) which is implemented as part of the [DEXSeq](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#3_Counting_reads) package. This is associated with differential exon usage analysis using DEXSeq and is activated when the parameter `--dexseq_exon` is enabled. Using the `--aggregation` parameter the pipeline will combine overlapping genes into a single aggregate gene. This approach can alternatively be skipped and any exons that overlap other exons from different genes will be skipped. Other important options to take note of are the `--alignment_quality` parameter which can be set by the user and defines the minimum alignment quality required for reads to be included (defined in 5th column of a given SAM file) (default: 10). Prior to quantification with HTSeq DEXSeq provides an annotation preparation script which takes a GTF file as input and returns a GFF file. Users may instead wish to define their own GFF file and skip this annotation preparation skip by supplying it using the `--gff_dexseq` parameter.

## Reference genome files

Like [nf-core/rnaseq](https://github.com/nf-core/rnaseq) the minimum reference genome requirements are a FASTA and GTF file. All other files required to run the pipeline can be generated from these files. If you wish to build new indices and save the output then the `--save_reference` parameter is required. It should be noted, however, that the sequence files and indices of many common species are available for download from [AWS iGenomes](https://nf-co.re/usage/reference_genomes)). Local copies of reference files may also be provided on the command line or via config files (e.g. `--salmon_index '/path/to/salmon_index.tar.gz'`). Where reference files are compressed (e.g. standard files with the `.gz` extension and indices folders with the `tar.gz` extension) these will be automatically uncompressed prior to use.

For nf-core/rnasplice to run it requires a FASTA file and GTF file. These can be specified manually, however users may wish to specify the more simple [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) approach:

- Using the `--genome` parameter followed by the genome build (e.g. GRCh37) will mean the FASTA and GTF file will be pulled from AWS-iGenomes, along with available indices required for a given analysis. Furthermore, a local download of AWS genomes from a previous run can also be used by changing the igenomes path with the `--igenomes_base` parameter.

If a GTF is not available a GFF may be used by specifying the `--gff` parameter. This will covert the GFF file into a GTF.

As in [nf-core/rnaseq](https://github.com/nf-core/rnaseq) if you are using a genome downloaded from AWS iGenomes and using `--aligner star_salmon` (default) the version of STAR to use for the alignment will be auto-detected (see [#808](https://github.com/nf-core/rnaseq/issues/808)).

Please note if you are using [GENCODE](https://www.gencodegenes.org/) reference genome files please specify the `--gencode` parameter. This is because reference files which come from GENCODE are different to ENSEMBL reference files and this can impact the running of the pipeline. Specifying this parameter can help to mitigate these differences. Furthermore it should be noted that when using GENCODE reference files if you are running Salmon, the `--gencode` flag will also be passed to the index building step (see [this issue](https://github.com/COMBINE-lab/salmon/issues/15)).

## Differential Exon Usage

Differential exon usage (**DEU**) can be completed by two different branches of the pipeline:

### DEXSeq DEU

Following HTSeq quantification you can estimate the differential exon usage using [DEXSeq](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html). This can be chosen by specifying the `--aligner star` or `--aligner star_salmon` followed by the `--dexseq_exon` parameter. Other important options to consider at this point for **DEU** with DEXSeq is the `--deu_lfc_denominator` parameter which defines the denominator of interest which DEXSeq will use to produce log2FoldChange calculations internally - the name of this denominator should match the name of one of the conditions from the samplesheet provided. If a user has multiple conditions a log2FoldChange will be calculated for each condition against this denominator.

### edgeR

Following featureCounts quantification differential exon usage can also be completed with [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).This can be chosen by specifying the `--aligner star` or `--aligner star_salmon` followed by the `--edger_exon` parameter. This DEU module will produce results for all possible comparisons of each condition against one another and you do not need to specify a denominator.

This module also produces differential exon expression results for each comparison.

## Differential Transcript Usage

### DEXSeq DTU

Differential transcript Usage (**DTU**) has been implemented from the workflow by [Love et al., 2018](https://f1000research.com/articles/7-952). This analysis can be chosen by specifying the `--aligner star_salmon` or `--pseudo_aligner salmon` followed by the `--dexseq_dtu` parameter. The parameter `--dtu_txi` should either be set as `--dtu_txi scaledTPM` or `--dtu_txi dtuScaledTPM` for successful **DTU** analysis.

Other important options to consider at this point for **DTU** with DEXSeq is the `--dtu_lfc_denominator` parameter which defines the denominator of interest which DEXSeq will use to produce log2FoldChange calculations internally - the name of this denominator should match the name of one of the conditions from the samplesheet provided.

Prior to DEXSeq DTU analysis filtering of genes and features with low expression is completed using [DRIMSeq](https://rdrr.io/bioc/DRIMSeq/man/dmFilter.html) which comes with a number of parameters which should be set by the user. By default these are all set to 0. For example, `--min_samps_gene_expr` defines the minimal number of samples where genes are required to be expressed, and `--min_gene_expr` defines the minimal level of gene expression for genes to be included in the downstream analysis. Similarly, `--min_samps_feature_expr` and `--min_samps_feature_prop` defines the minimal number of samples where features are required to be expressed at a minimal expression or proportion. This minimum level is further defined by additional filtering parameters `--min_feature_expr` and `--min_feature_prop` respectively. Further details of this filtering process can be see within the **DTU** workflow [here](https://f1000research.com/articles/7-952).

## Event based approaches

Two predominant event-based approaches have been implemented ([rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md) and [SUPPA](https://github.com/comprna/SUPPA)) in this pipeline and can be accessed through `--aligner` or `--pseudo_aligner` options:

### rMATS

[rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md) (replicate multivariate analysis of transcript splicing) is designed for detection of differential alternative splicing from replicate RNA-Seq data. If `--aligner` is set to `--aligner star` or `--aligner star_salmon` then rMATS can be accessed with the `--rmats` parameter.

At current, however, there are some restrictions to running an rMATS analysis:

- Samples need to have the same strandedness, read type (single-end/paired-end) and the samplesheet must have only one condition or two unique conditions. The rnasplice pipeline will automatically detect if you have a single condition and rMATS will run in single condition mode, otherwise a standard comparison will be run.
- The `--rmats_paired_stats` can be set to `true` only if there are two conditions and should not be run in single condition mode.

Furthermore, `--rmats_read_len` has to be set by the user and if the read length is variable, an average or median read length has to be specified.

### SUPPA2

You can run [SUPPA](https://github.com/comprna/SUPPA) for analyzing the splicing events across conditions following pseudo alignment and quantification with Salmon `--pseudo_aligner salmon` or after STAR alignment and Salmon quantification `--aligner star_salmon`, and when the `--suppa` parameter is supplied.

There are two main options for running an analysis with SUPPA - `--suppa_per_local_event` and `--suppa_per_isoform` (the latter is a DTU approach). When `--suppa_per_local_event` is set to `true`, local AS events are calculated and analyzed. When `--suppa_per_isoform` is set to `true`, transcript isoform events are calculated and analyzed.

#### Event Calculation

Events are calculated from user specified annotation files (e.g. GTF files). The parameter `--pool_genes` should be specified when creating ioe/ioi from annotations that are not loci-based. It should be noted that SUPPA advises users to utilse Ensembl and Gencode annotations to reduce errors at this stage of the analysis.

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

`PSI` files and `TPM` files are split based on the condition specified in metadata. e.g., condition1.psi, condition2.psi, condition1.tpm, condition2.tpm.
SUPPA then reads the `PSI` for the events and the transcript expression values from multiple samples, grouped by condition, and the `ioe`/`ioi` file, to calculate the events that are differentially spliced between a pair of conditions.

#### Cluster Analysis

Using `dpsi` file and `psivec` file, events are clustered according to `PSI` values across conditions.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/rnasplice --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rnasplice
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnasplice releases page](https://github.com/nf-core/rnasplice/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

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
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnasplice pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASPLICE:RNASPLICE:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASPLICE:RNASPLICE:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASPLICE:RNASPLICE:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

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
    withName: 'NFCORE_RNASPLICE:RNASPLICE:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASPLICE:RNASPLICE:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

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

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

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
