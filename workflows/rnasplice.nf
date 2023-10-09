/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['star', 'star_salmon'],
    pseudoaligners : ['salmon']
]

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRnasplice.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.fasta,
    params.gff,
    params.gtf,
    params.input,
    params.multiqc_config,
    params.salmon_index,
    params.star_index,
    params.transcript_fasta
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.contrasts) { ch_contrasts = file(params.contrasts) } else { exit 1, 'Input contrastsheet not specified!' }

// Check alignment parameters
def prepare_tool_indices  = []
if (!params.skip_alignment) { prepare_tool_indices << params.aligner        }
if (params.pseudo_aligner)  { prepare_tool_indices << params.pseudo_aligner }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { BEDTOOLS_GENOMECOV      } from '../modules/local/bedtools_genomecov'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                          } from '../subworkflows/local/input_check'
include { CONTRASTS_CHECK                                      } from '../subworkflows/local/contrasts_check'
include { PREPARE_GENOME                                       } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR                                           } from '../subworkflows/local/align_star'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_SALMON          } from '../subworkflows/local/tx2gene_tximport'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_STAR_SALMON     } from '../subworkflows/local/tx2gene_tximport'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_SALMON      } from '../subworkflows/local/drimseq_dexseq_dtu'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_STAR_SALMON } from '../subworkflows/local/drimseq_dexseq_dtu'
include { RMATS                                                } from '../subworkflows/local/rmats'
include { DEXSEQ_DEU                                           } from '../subworkflows/local/dexseq_deu'
include { EDGER_DEU                                            } from '../subworkflows/local/edger_deu'
include { SUPPA as SUPPA_SALMON                                } from '../subworkflows/local/suppa'
include { SUPPA as SUPPA_STAR_SALMON                           } from '../subworkflows/local/suppa'
include { VISUALISE_MISO                                       } from '../subworkflows/local/visualise_miso'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SALMON_QUANT as SALMON_QUANT_SALMON } from '../modules/nf-core/salmon/quant/main'
include { SALMON_QUANT as SALMON_QUANT_STAR   } from '../modules/nf-core/salmon/quant/main'
include { MULTIQC                             } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                           } from '../modules/nf-core/cat/fastq/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE                                } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BAM_SORT_STATS_SAMTOOLS                                         } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []
def pass_trimmed_reads = [:]

workflow RNASPLICE {

    // Create channel for software versions (will be added to throughout pipeline)
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.transcript_fasta,
        params.star_index,
        params.salmon_index,
        params.gff_dexseq,
        params.suppa_tpm,
        is_aws_igenome,
        prepare_tool_indices,
        params.source,
        params.gencode
    )
    ch_meta_fasta = PREPARE_GENOME.out.fasta.map { [ [:], it ] }
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    switch (params.source) {
        case 'fastq':
            INPUT_CHECK (
                ch_input,
                params.source
            )
            .reads
            .map {
                meta, fastq ->
                    new_id = meta.id - ~/_T\d+/
                    [ meta + [id: new_id], fastq ]
            }
            .groupTuple()
            .branch {
                meta, fastq ->
                    single  : fastq.size() == 1
                        return [ meta, fastq.flatten() ]
                    multiple: fastq.size() > 1
                        return [ meta, fastq.flatten() ]
            }
            .set { ch_fastq }
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
            break;
        case 'genome_bam':
            INPUT_CHECK (
                ch_input,
                params.source
            )
            .reads
            .set { ch_genome_bam }
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
            break;
        case 'transcriptome_bam':
            INPUT_CHECK (
                ch_input,
                params.source
            )
            .reads
            .set { ch_transcriptome_bam }
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
            break;
        case 'salmon_results':
            INPUT_CHECK (
                ch_input,
                params.source
            )
            .reads
            .set { ch_salmon_results }
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
            break;
    }

    // Create samplesheet channel (after input check)
    ch_samplesheet = Channel.fromPath(params.input)

    //
    // SUBWORKFLOW: Read in contrastsheet, validate and stage input files
    //
    CONTRASTS_CHECK (
        ch_contrasts
    )
    ch_versions = ch_versions.mix(CONTRASTS_CHECK.out.versions)

    // Create contrastsheet channel (after contrasts check)
    ch_contrastsheet = Channel.fromPath(params.contrasts)


    // Check rMATS parameters specified correctly
    if (params.rmats && params.source == 'fastq') {
            WorkflowRnasplice.rmatsReadError(INPUT_CHECK.out.reads, log)
            WorkflowRnasplice.rmatsStrandednessError(INPUT_CHECK.out.reads, log)
    }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    if (params.source == 'fastq') {
        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC and trim adapters with TrimGalore!
    //
    if (params.source == 'fastq') {

        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_cat_fastq,
            params.skip_fastqc,
            false,                    // run UMI false
            true,                     // skip UMI true
            params.skip_trimming,
            0,                        // No UMI discard
            0                         // No filtering 0 min_trimmed_reads
        )

        ch_trim_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_versions        = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)

        //
        // Get list of samples that failed trimming threshold for MultiQC report
        //
        ch_trim_read_count
            .map {
                meta, num_reads ->
                    pass_trimmed_reads[meta.id] = true
                    if (num_reads <= params.min_trimmed_reads.toFloat()) {
                        pass_trimmed_reads[meta.id] = false
                        return [ "$meta.id\t$num_reads" ]
                    }
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "Reads after trimming"]
                    WorkflowRnasplice.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_trimming_multiqc }

    }

    //
    // SUBWORKFLOW: Alignment with STAR
    //

    if (params.source == 'genome_bam') {

        BAM_SORT_STATS_SAMTOOLS ( ch_genome_bam, ch_meta_fasta )

        ch_genome_bam        = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_genome_bam_index  = BAM_SORT_STATS_SAMTOOLS.out.bai
        ch_samtools_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_samtools_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_samtools_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
        ch_versions          = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    }

    if (params.source == 'transcriptome_bam') {

        BAM_SORT_STATS_SAMTOOLS ( ch_transcriptome_bam, ch_meta_fasta )

        ch_transcriptome_bam        = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_transcriptome_bam_index  = BAM_SORT_STATS_SAMTOOLS.out.bai
        ch_samtools_stats           = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_samtools_flagstat        = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_samtools_idxstats        = BAM_SORT_STATS_SAMTOOLS.out.idxstats
        ch_versions                 = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    }

    if ((params.source == 'fastq') && !params.skip_alignment && ( params.aligner == 'star' || params.aligner == 'star_salmon')) {

        ALIGN_STAR (
            ch_trim_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            ch_meta_fasta
        )

        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final

        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }

        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    if ((params.source == 'genome_bam') || (params.source == 'fastq') && (!params.skip_alignment && ( params.aligner == 'star' || params.aligner == 'star_salmon'))) {

        //
        // SUBWORKFLOW: Run DEXSeq DEU
        //

        if (params.dexseq_exon) {

            ch_dexseq_gff = params.gff_dexseq ? PREPARE_GENOME.out.dexseq_gff : ""

            DEXSEQ_DEU(
                PREPARE_GENOME.out.gtf,
                ch_genome_bam,
                ch_dexseq_gff,
                ch_samplesheet,
                ch_contrastsheet,
                params.n_dexseq_plot,
                params.aggregation,
                params.alignment_quality
            )

            ch_versions = ch_versions.mix(DEXSEQ_DEU.out.versions)
        }

        //
        // SUBWORKFLOW: Run edgeR DEU
        //

        if (params.edger_exon) {

            EDGER_DEU(
                PREPARE_GENOME.out.gtf,
                ch_genome_bam,
                ch_samplesheet,
                ch_contrastsheet,
                params.n_edger_plot
            )

            ch_versions = ch_versions.mix(EDGER_DEU.out.versions)

        }

        //
        // SUBWORKFLOW: Run rMATS
        //

        if (params.rmats) {

            //
            // Create channel grouped by condition: [ [condition1, [condition1_metas], [group1_bams]], [condition2, [condition2_metas], [condition2_bams]]] if there are more than one condition
            //

            // If genome bam provided, mimic channel created by ALIGN_STAR
            if (params.source == 'genome_bam') {

                BAM_SORT_STATS_SAMTOOLS
                        .out
                        .bam
                        .map { meta, bam -> [meta.condition, meta, bam] }
                        .groupTuple(by:0)
                        .set { ch_genome_bam_conditions }

            } else {

                ALIGN_STAR
                    .out
                    .bam
                    .map { meta, bam -> [meta.condition, meta, bam] }
                    .groupTuple(by:0)
                    .set { ch_genome_bam_conditions }

            }

            //
            // Create variable to check if samples have one condition or two
            //

            is_single_condition = WorkflowRnasplice.isSingleCondition(ch_input)

            //
            // SUBWORKFLOW: Run rMATS
            //

            RMATS (
                ch_contrastsheet,
                ch_genome_bam_conditions,
                PREPARE_GENOME.out.gtf,
                is_single_condition,
                params.rmats_read_len,
                params.rmats_splice_diff_cutoff,
                params.rmats_novel_splice_site,
                params.rmats_min_intron_len,
                params.rmats_max_exon_len,
                params.rmats_paired_stats
            )

            ch_versions = ch_versions.mix(RMATS.out.versions)

        }

        if (params.sashimi_plot == true) {

            VISUALISE_MISO (
                PREPARE_GENOME.out.gtf,
                ch_genome_bam,
                ch_genome_bam_index,
                params.miso_read_len,
                params.fig_width,
                params.fig_height,
                params.miso_genes,
                params.miso_genes_file ?: false,
            )

            ch_versions = ch_versions.mix(VISUALISE_MISO.out.versions)

        }

    }

    //
    // SUBWORKFLOW: Count reads from BAM alignments using Salmon
    //

    // If needed run Salmon on transcriptome bam

    if ((params.source == 'transcriptome_bam') || (params.source == 'fastq') && (!params.skip_alignment && params.aligner == 'star_salmon')) {

        alignment_mode = true
        ch_salmon_index = ch_dummy_file

        // Run Salmon quant, Run tx2gene.py (tx2gene for Salmon txImport Quantification), then finally runs tximport
        SALMON_QUANT_STAR (
            ch_transcriptome_bam,
            ch_salmon_index,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.transcript_fasta,
            alignment_mode,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(SALMON_QUANT_STAR.out.versions)

        // Collect Salmon quant output
        ch_salmon_results = SALMON_QUANT_STAR.out.results

        //
        // SUBWORKFLOW: Run Tximport and produce tx2gene from gtf using gffread
        //
        TX2GENE_TXIMPORT_STAR_SALMON (
            ch_salmon_results,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(TX2GENE_TXIMPORT_STAR_SALMON.out.versions)

        //
        // SUBWORKFLOW: Run Dexseq DTU
        //

        if (params.dexseq_dtu) {

            if (params.dtu_txi == "dtuScaledTPM") {

                ch_txi = TX2GENE_TXIMPORT_STAR_SALMON.out.txi_dtu

            } else if (params.dtu_txi == "scaledTPM") {

                ch_txi = TX2GENE_TXIMPORT_STAR_SALMON.out.txi_s

            }

            DRIMSEQ_DEXSEQ_DTU_STAR_SALMON (
                ch_txi,
                TX2GENE_TXIMPORT_STAR_SALMON.out.tximport_tx2gene,
                ch_samplesheet,
                ch_contrastsheet,
                params.n_dexseq_plot,
                params.min_samps_gene_expr,
                params.min_samps_feature_expr,
                params.min_samps_feature_prop,
                params.min_feature_expr,
                params.min_feature_prop,
                params.min_gene_expr
            )

            ch_versions = ch_versions.mix(DRIMSEQ_DEXSEQ_DTU_STAR_SALMON.out.versions)

        }

        //
        // SUBWORKFLOW: SUPPA
        //

        if (params.suppa) {

            // Get Suppa tpm either from tximport or user supplied
            ch_suppa_tpm = params.suppa_tpm ? PREPARE_GENOME.out.suppa_tpm : TX2GENE_TXIMPORT_STAR_SALMON.out.suppa_tpm

            // Run SUPPA
            SUPPA_STAR_SALMON (
                PREPARE_GENOME.out.gtf,
                ch_suppa_tpm,
                ch_samplesheet,
                ch_contrastsheet,
                params.suppa_per_local_event,
                params.generateevents_boundary,
                params.generateevents_threshold,
                params.generateevents_exon_length,
                params.generateevents_event_type,
                params.generateevents_pool_genes,
                params.psiperevent_total_filter,
                params.diffsplice_local_event,
                params.diffsplice_method,
                params.diffsplice_area,
                params.diffsplice_lower_bound,
                params.diffsplice_alpha,
                params.diffsplice_tpm_threshold,
                params.diffsplice_nan_threshold,
                params.diffsplice_gene_correction,
                params.diffsplice_paired,
                params.diffsplice_median,
                params.clusterevents_local_event,
                params.clusterevents_dpsithreshold,
                params.clusterevents_eps,
                params.clusterevents_metric,
                params.clusterevents_min_pts,
                params.clusterevents_method,
                params.clusterevents_sigthreshold ?: false,
                params.clusterevents_separation ?: false,
                params.suppa_per_isoform
            )

            ch_versions = ch_versions.mix(SUPPA_STAR_SALMON.out.versions)

        }

    }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //

    if (params.source == 'fastq' && params.pseudo_aligner == 'salmon') {

        alignment_mode = false
        ch_transcript_fasta = ch_dummy_file

        SALMON_QUANT_SALMON (
            ch_trim_reads,
            PREPARE_GENOME.out.salmon_index,
            PREPARE_GENOME.out.gtf,
            ch_transcript_fasta,
            alignment_mode,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(SALMON_QUANT_SALMON.out.versions)

        // Collect Salmon quant output
        ch_salmon_results = SALMON_QUANT_SALMON.out.results

        //
        // SUBWORKFLOW: Run Tximport and produce tx2gene from gtf using gffread
        //

        TX2GENE_TXIMPORT_SALMON (
            ch_salmon_results,
            PREPARE_GENOME.out.gtf
        )

        ch_versions = ch_versions.mix(TX2GENE_TXIMPORT_SALMON.out.versions)


    } else if ( params.source == 'salmon_results') {

        TX2GENE_TXIMPORT_SALMON (
            ch_salmon_results,
            PREPARE_GENOME.out.gtf
        )

        ch_versions = ch_versions.mix(TX2GENE_TXIMPORT_SALMON.out.versions)

    }


    if ((params.pseudo_aligner == 'salmon' && params.source == 'fastq') || (params.source == 'salmon_results')) {

        //
        // SUBWORKFLOW: Run Dexseq DTU
        //
        if (params.dexseq_dtu) {
            if (params.dtu_txi == "dtuScaledTPM") {
                ch_txi = TX2GENE_TXIMPORT_SALMON.out.txi_dtu
            } else if (params.dtu_txi == "scaledTPM") {
                ch_txi = TX2GENE_TXIMPORT_SALMON.out.txi_s
            }

            DRIMSEQ_DEXSEQ_DTU_SALMON (
                ch_txi,
                TX2GENE_TXIMPORT_SALMON.out.tximport_tx2gene,
                ch_samplesheet,
                ch_contrastsheet,
                params.n_dexseq_plot,
                params.min_samps_gene_expr,
                params.min_samps_feature_expr,
                params.min_samps_feature_prop,
                params.min_feature_expr,
                params.min_feature_prop,
                params.min_gene_expr
            )
            ch_versions = ch_versions.mix(DRIMSEQ_DEXSEQ_DTU_SALMON.out.versions)
        }

        //
        // SUBWORKFLOW: SUPPA
        //
        if (params.suppa) {

            // Get Suppa tpm either from tximport or user supplied
            ch_suppa_tpm = params.suppa_tpm ? PREPARE_GENOME.out.suppa_tpm : TX2GENE_TXIMPORT_SALMON.out.suppa_tpm

            // Run SUPPA
            SUPPA_SALMON (
                PREPARE_GENOME.out.gtf,
                ch_suppa_tpm,
                ch_samplesheet,
                ch_contrastsheet,
                params.suppa_per_local_event,
                params.generateevents_boundary,
                params.generateevents_threshold,
                params.generateevents_exon_length,
                params.generateevents_event_type,
                params.generateevents_pool_genes,
                params.psiperevent_total_filter,
                params.diffsplice_local_event,
                params.diffsplice_method,
                params.diffsplice_area,
                params.diffsplice_lower_bound,
                params.diffsplice_alpha,
                params.diffsplice_tpm_threshold,
                params.diffsplice_nan_threshold,
                params.diffsplice_gene_correction,
                params.diffsplice_paired,
                params.diffsplice_median,
                params.clusterevents_local_event,
                params.clusterevents_dpsithreshold,
                params.clusterevents_eps,
                params.clusterevents_metric,
                params.clusterevents_min_pts,
                params.clusterevents_method,
                params.clusterevents_sigthreshold ?: false,
                params.clusterevents_separation ?: false,
                params.suppa_per_isoform
            )

            ch_versions = ch_versions.mix(SUPPA_SALMON.out.versions)

        }

    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (((params.source == 'genome_bam') && !params.skip_bigwig) || (!params.skip_alignment && !params.skip_bigwig && params.source == 'fastq')) {

        BEDTOOLS_GENOMECOV (
            ch_genome_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

	    //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_TO_BIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV.out.bedgraph_forward,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDGRAPH_TO_BIGWIG_FORWARD.out.versions)

        BEDGRAPH_TO_BIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV.out.bedgraph_reverse,
            PREPARE_GENOME.out.chrom_sizes
        )
    }


    //
    // MODULE: Collect version information across pipeline
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowRnasplice.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRnasplice.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    if (params.source == 'fastq') {

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]))

    }

    if (params.pseudo_aligner == 'salmon' && params.source == 'fastq'){

        ch_multiqc_files = ch_multiqc_files.mix(ch_salmon_results.collect{it[1]}.ifEmpty([]))

    }

    if (!params.skip_alignment && params.source == 'fastq' && (params.aligner == 'star_salmon' || params.aligner == "star")){

        ch_multiqc_files = ch_multiqc_files.mix(ch_star_multiqc.collect{it[1]}.ifEmpty([]))

    }

    if ((params.source == 'genome_bam' || params.source == 'transcriptome_bam') || (!params.skip_alignment && params.source == 'fastq' && (params.aligner == 'star_salmon' || params.aligner == "star"))){

        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))


        if (params.edger_exon && params.source != 'transcriptome_bam') {

            ch_multiqc_files = ch_multiqc_files.mix(EDGER_DEU.out.featureCounts_summary.collect{it[1]}.ifEmpty([]))

        }
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
