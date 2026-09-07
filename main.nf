#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnasplice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnasplice
    Website: https://nf-co.re/rnasplice
    Slack  : https://nfcore.slack.com/channels/rnasplice
----------------------------------------------------------------------------------------
*/

params {

    // Input options
    input: String
    contrasts: String
    source: String = 'fastq'

    // References
    fasta: String = getGenomeAttribute('fasta')
    gtf: String? = getGenomeAttribute('gtf')
    gff: String? = getGenomeAttribute('gff')
    transcript_fasta: String
    gtf_extra_attributes: String = 'gene_name'
    gtf_group_features: String = 'gene_id'
    gencode: Boolean
    save_reference: Boolean
    igenomes_base: String = 's3://ngi-igenomes/igenomes/'

    // Trimming
    clip_r1: Integer
    clip_r2: Integer
    three_prime_clip_r1: Integer
    three_prime_clip_r2: Integer
    trim_nextseq: Integer
    save_trimmed: Boolean
    skip_trimming: Boolean
    skip_trimgalore_fastqc: Boolean
    min_trimmed_reads: Integer = 10000

    // Alignment
    aligner: String = 'star_salmon'
    pseudo_aligner: String = 'salmon'
    bam_csi_index: Boolean
    seq_center: String
    salmon_quant_libtype: String
    star_index: String? = getGenomeAttribute('star')
    salmon_index: String? = getGenomeAttribute('salmon')
    star_ignore_sjdbgtf: Boolean
    skip_alignment: Boolean
    save_unaligned: Boolean
    save_align_intermeds: Boolean
    save_merged_fastq: Boolean

    // QC
    skip_bigwig: Boolean = true
    skip_fastqc: Boolean

    // rMATs
    rmats: Boolean = true
    rmats_splice_diff_cutoff: Float = 0.0001
    rmats_paired_stats: Boolean = true
    rmats_read_len: Integer = 40
    rmats_novel_splice_site: Boolean
    rmats_min_intron_len: Integer = 50
    rmats_max_exon_len: Integer = 500

    // DEXSeq DEU
    dexseq_exon: Boolean = true
    save_dexseq_annotation: Boolean
    gff_dexseq: String
    alignment_quality: Integer = 10
    aggregation: Boolean = true
    save_dexseq_plot: Boolean = true
    n_dexseq_plot: Integer = 10
    ignore_tx_version: Boolean = true

    // edgeR DEU
    edger_exon: Boolean = true
    save_edger_plot: Boolean = true
    n_edger_plot: Integer = 10

    // DEXSeq DTU
    dexseq_dtu: Boolean = true
    dtu_txi: String = 'dtuScaledTPM'

    // Miso
    sashimi_plot: Boolean = true
    miso_genes: String = 'ENSG00000004961, ENSG00000005302, ENSG00000147403'
    miso_genes_file: String
    miso_read_len: Integer = 75
    fig_height: Integer = 7
    fig_width: Integer = 7

    // Leafcutter
    leafcutter: Boolean

    // DRIMSeq Filtering
    min_samps_feature_expr: Integer = 2
    min_samps_feature_prop: Integer = 2
    min_samps_gene_expr: Integer = 4
    min_gene_expr: Integer = 10
    min_feature_expr: Integer = 10
    min_feature_prop: Float = 0.1

    // SUPPA options
    suppa: Boolean = true
    suppa_per_local_event: Boolean = true
    suppa_per_isoform: Boolean = true
    suppa_tpm: String

    // SUPPA generateEvents options
    generateevents_pool_genes: Boolean = true
    generateevents_event_type: String = 'SE SS MX RI FL'
    generateevents_boundary: String = 'S'
    generateevents_threshold: Integer = 10
    generateevents_exon_length: Integer = 100
    psiperevent_total_filter: Integer = 0

    // SUPPA Diffsplice options
    diffsplice_local_event: Boolean = true
    diffsplice_isoform: Boolean = true
    diffsplice_method: String = 'empirical'
    diffsplice_area: Integer = 1000
    diffsplice_lower_bound: Integer = 0
    diffsplice_gene_correction: Boolean = true
    diffsplice_paired: Boolean = true
    diffsplice_alpha: Float = 0.05
    diffsplice_median: Boolean
    diffsplice_tpm_threshold: Integer = 0
    diffsplice_nan_threshold: Integer = 0

    // SUPPA Cluster options
    clusterevents_local_event: Boolean = true
    clusterevents_isoform: Boolean = true
    clusterevents_sigthreshold: Float
    clusterevents_dpsithreshold: Float = 0.05
    clusterevents_eps: Float = 0.05
    clusterevents_metric: String = 'euclidean'
    clusterevents_separation: Float
    clusterevents_min_pts: Integer = 20
    clusterevents_method: String = 'DBSCAN'

    // IsoformSwitchAnalyzer options
    isoformswitchanalyzer: Boolean = true
    isoformswitchanalyzer_alpha: Float = 0.05
    isoformswitchanalyzer_dIF: Float = 0.1

    // MultiQC options
    multiqc_config: String
    multiqc_logo: String
    multiqc_title: String
    max_multiqc_email_size: String = 25.MB
    multiqc_methods_description: String

    // Boilerplate options
    outdir: String
    publish_dir_mode: String = 'copy'
    email: String
    email_on_fail: String
    plaintext_email: Boolean
    help: Boolean
    help_full: Boolean
    show_hidden: Boolean
    version: Boolean

    // Schema validation default options
    validate_params: Boolean = true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASPLICE               } from './workflows/rnasplice'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_RNASPLICE {

    main:

    def is_aws_igenome = params.fasta && params.gtf && (file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')

    //
    // SUBWORKFLOW: Prepare reference genome files
    //

    PREPARE_GENOME(
        params.fasta,
        params.gtf,
        params.gff,
        params.transcript_fasta,
        params.star_index,
        params.salmon_index,
        params.gff_dexseq,
        params.suppa_tpm,
        params.gencode,
        is_aws_igenome,
    )

    ch_samplesheet = channel.value(file(params.input, checkIfExists: true))
    ch_contrastsheet = channel.value(file(params.contrasts, checkIfExists: true))

    //
    // WORKFLOW: Run pipeline
    //
    RNASPLICE(
        ch_samplesheet,
        ch_contrastsheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.dexseq_gff,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.suppa_tpm,
        PREPARE_GENOME.out.chrom_sizes,
        is_aws_igenome,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    multiqc_report = RNASPLICE.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNASPLICE()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_RNASPLICE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}
