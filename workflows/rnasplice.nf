/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['star', 'star_salmon'],
    pseudoaligners : ['salmon']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnasplice.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta,
    params.gtf, params.gff,
    params.star_index, params.salmon_index
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepare_tool_indices  = []
if (!params.skip_alignment) { prepare_tool_indices << params.aligner        }
if (params.pseudo_aligner)  { prepare_tool_indices << params.pseudo_aligner }

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false

if (params.fasta && params.gtf) {

    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}
// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

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
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore'
include { ALIGN_STAR        } from '../subworkflows/local/align_star'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_SALMON      } from '../subworkflows/local/tx2gene_tximport'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_STAR_SALMON } from '../subworkflows/local/tx2gene_tximport'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_SALMON } from '../subworkflows/local/drimseq_dexseq_dtu'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_STAR_SALMON } from '../subworkflows/local/drimseq_dexseq_dtu'
include { RMATS             } from '../subworkflows/local/rmats'
include { DEXSEQ_DEU        } from '../subworkflows/local/dexseq_deu'
include { EDGER_DEU         } from '../subworkflows/local/edger_deu'
include { SUPPA as SUPPA_SALMON       } from '../subworkflows/local/suppa'
include { SUPPA as SUPPA_STAR_SALMON  } from '../subworkflows/local/suppa'

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
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD       } from '../subworkflows/nf-core/bedgraph_to_bigwig'
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE       } from '../subworkflows/nf-core/bedgraph_to_bigwig'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

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
        prepare_tool_indices
    )

    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    // Run Input check subworkflow
    INPUT_CHECK (
	    ch_input
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    // Create samplesheet channel (after input check)
    ch_samplesheet = Channel.fromPath(params.input)

    // Take software versions from input check (.first() not required)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Check rMATS parameters specified correctly
    if (params.rmats) {

        WorkflowRnasplice.rmatsReadError(INPUT_CHECK.out.reads, log)

        WorkflowRnasplice.rmatsStrandednessError(INPUT_CHECK.out.reads, log)

        WorkflowRnasplice.rmatsConditionError(INPUT_CHECK.out.reads, log)

    }

    // Check DEXSeq parameters specified correctly

    if (params.dexseq_exon || params.dexseq_dtu) {

        WorkflowRnasplice.denominatorExistsError(params, log, ch_input)

    }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Read QC, and trimming
    //

    // Run FastQC and TrimGalore
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )

    // Take software versions from subworkflow (.first() not required)
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    // Collect trimmed reads from Trimgalore
    ch_trim_reads = FASTQC_TRIMGALORE.out.reads

    //
    // SUBWORKFLOW: Alignment with STAR
    //

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()

    if (!params.skip_alignment && ( params.aligner == 'star' || params.aligner == 'star_salmon')) {

        ALIGN_STAR (
            ch_trim_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome
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

        //
        // SUBWORKFLOW: Run DEXSeq DEU branch (params.dexseq_exon = true)
        //

        if (params.dexseq_exon) {

            ch_dexseq_gff = params.gff_dexseq ? PREPARE_GENOME.out.dexseq_gff : ""
            def read_method = "htseq"

            DEXSEQ_DEU(
                PREPARE_GENOME.out.gtf,
                ch_genome_bam,
                ch_dexseq_gff,
                ch_samplesheet,
                read_method
            )

            ch_versions = ch_versions.mix(DEXSEQ_DEU.out.versions)
        }

        //
        // SUBWORKFLOW: Run edgeR DEU branch (params.edger_exon = true)
        //

        if (params.edger_exon) {

            EDGER_DEU(
                PREPARE_GENOME.out.gtf,
                ch_genome_bam,
                ch_samplesheet
            )

            ch_versions = ch_versions.mix(EDGER_DEU.out.versions)

        }

        //
        // Run rMATS subworkflow if rmats parameter true:
        //

        if (params.rmats) {

            //
            // Create channel grouped by condition: [ [condition1, [condition1_metas], [group1_bams]], [condition2, [condition2_metas], [condition2_bams]]] if there are more than one condition
            //

            ALIGN_STAR
                .out
                .bam
                .map { meta, bam -> [meta.condition, meta, bam] }
                .groupTuple(by:0)
                .toSortedList()
                .set { ch_genome_bam_conditions }

            //
            // Create variable to check if samples have one condition or two
            //

            is_single_condition = WorkflowRnasplice.isSingleCondition(ch_input)

            //
            // SUBWORKFLOW: Run rMATS
            //

            RMATS (
                ch_genome_bam_conditions,
                PREPARE_GENOME.out.gtf,
                is_single_condition
            )

            ch_versions = ch_versions.mix(RMATS.out.versions)

        }

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //

        // If needed run Salmon on transcriptome bam
        if (!params.skip_alignment && params.aligner == 'star_salmon') {

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


            //
            // SUBWORKFLOW: Run Tximport and produce tx2gene from gtf using gffread
            //

            TX2GENE_TXIMPORT_STAR_SALMON (
                SALMON_QUANT_STAR.out.results.collect{it[1]},
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
                    ch_samplesheet
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
                )

                ch_versions = ch_versions.mix(SUPPA_STAR_SALMON.out.versions)

            }


        }

    }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //

    if (params.pseudo_aligner == 'salmon') {

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

        // Collect Salmon quant output
        ch_salmon_multiqc = SALMON_QUANT_SALMON.out.results

        // Take software versions from subworkflow (.first() not required)
        ch_versions = ch_versions.mix(SALMON_QUANT_SALMON.out.versions)

        //
        // SUBWORKFLOW: Run Tximport and produce tx2gene from gtf using gffread
        //

        TX2GENE_TXIMPORT_SALMON (
            SALMON_QUANT_SALMON.out.results.collect{it[1]},
            PREPARE_GENOME.out.gtf
        )

        ch_versions = ch_versions.mix(TX2GENE_TXIMPORT_SALMON.out.versions)

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
                ch_samplesheet
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
            )

            ch_versions = ch_versions.mix(SUPPA_SALMON.out.versions)

        }
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {

        BEDTOOLS_GENOMECOV (
            ch_genome_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

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

    methods_description    = WorkflowRnasplice.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    if (params.pseudo_aligner == 'salmon'){

        ch_multiqc_files = ch_multiqc_files.mix(ch_salmon_multiqc.collect{it[1]}.ifEmpty([]))

    }

    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == "star")){

        ch_multiqc_files = ch_multiqc_files.mix(ch_star_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))

        if (params.edger_exon) {

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
