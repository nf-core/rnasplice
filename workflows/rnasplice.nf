/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnasplice.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }
if (params.pseudo_aligner)  { prepareToolIndices << params.pseudo_aligner }

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

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { BEDTOOLS_GENOMECOV      } from '../modules/local/bedtools_genomecov'
include { STAR_ALIGN_IGENOMES     } from '../modules/local/star_align_igenomes'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore'
include { TX2GENE_TXIMPORT as SALMON_TX2GENE_TXIMPORT      } from '../subworkflows/local/tx2gene_tximport'
include { TX2GENE_TXIMPORT as STAR_SALMON_TX2GENE_TXIMPORT } from '../subworkflows/local/tx2gene_tximport'
include { DRIMSEQ_DEXSEQ_DTU as SALMON_DEXSEQ_DTU } from '../subworkflows/local/drimseq_dexseq_dtu'
include { DRIMSEQ_DEXSEQ_DTU as STAR_SALMON_DEXSEQ_DTU } from '../subworkflows/local/drimseq_dexseq_dtu'
include { RMATS             } from '../subworkflows/local/rmats'
include { DEXSEQ_DEU        } from '../subworkflows/local/dexseq_deu'
include { EDGER_DEU         } from '../subworkflows/local/edger_deu'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SALMON_QUANT                      } from '../modules/nf-core/salmon/quant/main'
include { SALMON_QUANT as STAR_SALMON_QUANT } from '../modules/nf-core/salmon/quant/main'
include { STAR_ALIGN                        } from '../modules/nf-core/star/align/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                         } from '../modules/nf-core/cat/fastq/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { BAM_SORT_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_samtools'
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

    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type

    PREPARE_GENOME (
        prepareToolIndices,
        biotype,
        is_aws_igenome
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

    // Take software versions from input check (.first() not required)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    if (params.rmats) {

        INPUT_CHECK
            .out
            .reads
            .map { meta, fastq -> meta.single_end}
            .unique()
            .collect()
            .map {
                if(it.size() > 1) {
                    exit 1, "ERROR: Please check input samplesheet -> Cannot run rMats with mixed single and paired end samples"
                } else {
                    return it
                }
            }

        INPUT_CHECK
            .out
            .reads
            .map { meta, fastq -> meta.strandedness}
            .unique()
            .collect()
            .map {
                if(it.size() > 1) {
                    exit 1, "ERROR: Please check input samplesheet -> Cannot run rMats with mixed stranded samples"
                } else {
                    return it
                }
            }
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
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    
    // Take software versions from subworkflow (.first() not required)
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    // Collect trimmed reads from Trimgalore
    ch_trim_reads = FASTQC_TRIMGALORE.out.reads

    //
    // MODULE: Align reads using STAR
    //
            
    // Align with STAR
    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == "star")) {
        
        if (is_aws_igenome) {

            STAR_ALIGN_IGENOMES ( 
                ch_trim_reads, 
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.gtf, 
                params.star_ignore_sjdbgtf,
                params.seq_platform ?: '',
                params.seq_center ?: '' 
            )
            // Collect STAR_ALIGN_IGENOMES output 
            ch_orig_bam          = STAR_ALIGN_IGENOMES.out.bam            // channel: [ val(meta), bam            ]
            ch_log_final         = STAR_ALIGN_IGENOMES.out.log_final      // channel: [ val(meta), log_final      ]
            ch_log_out           = STAR_ALIGN_IGENOMES.out.log_out        // channel: [ val(meta), log_out        ]
            ch_log_progress      = STAR_ALIGN_IGENOMES.out.log_progress   // channel: [ val(meta), log_progress   ]
            ch_bam_sorted        = STAR_ALIGN_IGENOMES.out.bam_sorted     // channel: [ val(meta), bam_sorted     ]
            ch_transcriptome_bam = STAR_ALIGN_IGENOMES.out.bam_transcript // channel: [ val(meta), bam_transcript ]
            ch_fastq             = STAR_ALIGN_IGENOMES.out.fastq          // channel: [ val(meta), fastq          ]   
            ch_tab               = STAR_ALIGN_IGENOMES.out.tab            // channel: [ val(meta), tab            ]

            // Collect software version
            ch_versions       = ch_versions.mix(STAR_ALIGN_IGENOMES.out.versions.first())

        } else {

            // Run Star alignment module 
            // currently params.seq_platform not specified in the nextflow.config in rnaseq code base
            STAR_ALIGN ( 
                ch_trim_reads, 
                PREPARE_GENOME.out.star_index, 
                PREPARE_GENOME.out.gtf, 
                params.star_ignore_sjdbgtf, 
                params.seq_platform ?: '', 
                params.seq_center ?: ''
            )
            // Collect STAR output 
            ch_orig_bam          = STAR_ALIGN.out.bam            // channel: [ val(meta), bam            ]
            ch_log_final         = STAR_ALIGN.out.log_final      // channel: [ val(meta), log_final      ]
            ch_log_out           = STAR_ALIGN.out.log_out        // channel: [ val(meta), log_out        ]
            ch_log_progress      = STAR_ALIGN.out.log_progress   // channel: [ val(meta), log_progress   ]
            ch_bam_sorted        = STAR_ALIGN.out.bam_sorted     // channel: [ val(meta), bam_sorted     ]
            ch_transcriptome_bam = STAR_ALIGN.out.bam_transcript // channel: [ val(meta), bam_transcript ]
            ch_fastq             = STAR_ALIGN.out.fastq          // channel: [ val(meta), fastq          ]   
            ch_tab               = STAR_ALIGN.out.tab            // channel: [ val(meta), tab            ]

            // Collect software version
            ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())
        }

        //
        // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
        //

        // Run Samtools subworkflow (sort, index with stats)
        BAM_SORT_SAMTOOLS ( ch_orig_bam )

        // Collect Samtools output - sorted bam, indices (bai, csi)
        ch_genome_bam        = BAM_SORT_SAMTOOLS.out.bam         // channel: [ val(meta), [ bam ] ]
        ch_genome_bam_index  = BAM_SORT_SAMTOOLS.out.bai         // channel: [ val(meta), [ bai ] ]
        ch_genom_csi         = BAM_SORT_SAMTOOLS.out.csi         // channel: [ val(meta), [ csi ] ]

        // Collect Samtools stats output - stats, flagstat, idxstats
        ch_samtools_stats    = BAM_SORT_SAMTOOLS.out.stats       // channel: [ val(meta), [ stats ] ]
        ch_samtools_flagstat = BAM_SORT_SAMTOOLS.out.flagstat    // channel: [ val(meta), [ flagstat ] ]
        ch_samtools_idxstats = BAM_SORT_SAMTOOLS.out.idxstats    // channel: [ val(meta), [ idxstats ] ]
        
        // Collect software version
        ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

        //
        // SUBWORKFLOW: Run DEXSeq DEU branch (params.dexseq_exon = true)
        //

        if (params.dexseq_exon) {

            ch_dexseq_gff = params.gff_dexseq ? PREPARE_GENOME.out.dexseq_gff : ""
            ch_samplesheet = Channel.fromPath(params.input)
            def read_method = "htseq"
        
            DEXSEQ_DEU(
                PREPARE_GENOME.out.gtf,
                ch_genome_bam,
                ch_dexseq_gff,
                ch_samplesheet,
                read_method
            )
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

        }

        //
        // Run rMATS subworkflow if rmats paramater true:
        //

        if (params.rmats) {

            //
            // Create channel grouped by condition: [ [condition1, [condition1_metas], [group1_bams]], [condition2, [condition2_metas], [condition2_bams]]]
            //
                
            BAM_SORT_SAMTOOLS
                .out
                .bam
                .map { meta, bam -> [meta.condition, meta, bam] }
                .groupTuple(by:0)
                .toSortedList()
                .set { ch_genome_bam_conditions }
            
            //
            // SUBWORKFLOW: Run rMATS
            //
            
            RMATS ( 
                ch_genome_bam_conditions,
                PREPARE_GENOME.out.gtf 
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
            STAR_SALMON_QUANT (
                ch_transcriptome_bam,
                ch_salmon_index,
                PREPARE_GENOME.out.gtf,
                PREPARE_GENOME.out.transcript_fasta,
                alignment_mode,
                params.salmon_quant_libtype ?: ''
            )

            ch_versions = ch_versions.mix(STAR_SALMON_QUANT.out.versions)


            //
            // SUBWORKFLOW: Run Tximport and produce tx2gene from gtf using gffread
            //

            STAR_SALMON_TX2GENE_TXIMPORT (
                STAR_SALMON_QUANT.out.results.collect{it[1]},
                PREPARE_GENOME.out.gtf
            )

            //
            // SUBWORKFLOW: Run Dexseq DTU
            //
            
            if (params.dexseq_dtu) {

                ch_samplesheet = Channel.fromPath(params.input)

                if (params.dtu_txi == "dtuScaledTPM") {

                    ch_txi = SALMON_TX2GENE_TXIMPORT.out.txi_dtu

                } else if (params.dtu_txi == "scaledTPM") {

                    ch_txi = SALMON_TX2GENE_TXIMPORT.out.txi_s

                }

                STAR_SALMON_DEXSEQ_DTU (
                    ch_txi,
                    SALMON_TX2GENE_TXIMPORT.out.tximport_tx2gene,
                    ch_samplesheet
                )
            }
        }

    }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //

    if (params.pseudo_aligner == 'salmon') {
        
        alignment_mode = false
        ch_transcript_fasta = ch_dummy_file

        SALMON_QUANT (
            ch_trim_reads,
            PREPARE_GENOME.out.salmon_index,
            PREPARE_GENOME.out.gtf,
            ch_transcript_fasta,
            alignment_mode,
            params.salmon_quant_libtype ?: ''
        )

        // Collect Salmon quant output
        ch_salmon_multiqc = SALMON_QUANT.out.results

        // Take software versions from subworkflow (.first() not required)
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

        //
        // SUBWORKFLOW: Run Tximport and produce tx2gene from gtf using gffread
        //

        SALMON_TX2GENE_TXIMPORT (
            SALMON_QUANT.out.results.collect{it[1]},
            PREPARE_GENOME.out.gtf
        )

        //
        // SUBWORKFLOW: Run Dexseq DTU
        //
        if (params.dexseq_dtu) {

            ch_samplesheet = Channel.fromPath(params.input)

            if (params.dtu_txi == "dtuScaledTPM") {

                ch_txi = SALMON_TX2GENE_TXIMPORT.out.txi_dtu

            } else if (params.dtu_txi == "scaledTPM") {

                ch_txi = SALMON_TX2GENE_TXIMPORT.out.txi_s

            }

            SALMON_DEXSEQ_DTU (
                ch_txi,
                SALMON_TX2GENE_TXIMPORT.out.tximport_tx2gene,
                ch_samplesheet
            )
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

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    if (params.pseudo_aligner == 'salmon'){
        ch_multiqc_files = ch_multiqc_files.mix(ch_salmon_multiqc.collect{it[1]}.ifEmpty([]))
    }

    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == "star")){

        ch_multiqc_files = ch_multiqc_files.mix(ch_log_final.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))

    }

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
