/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FORWARD                } from '../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REVERSE                } from '../modules/nf-core/bedtools/genomecov'
include { ISOFORMSWITCHANALYZER                                           } from '../modules/local/isoformswitchanalyzer'
include { ALIGN_STAR                                                      } from '../subworkflows/local/align_star'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_SALMON                     } from '../subworkflows/local/tx2gene_tximport'
include { TX2GENE_TXIMPORT as TX2GENE_TXIMPORT_STAR_SALMON                } from '../subworkflows/local/tx2gene_tximport'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_SALMON                 } from '../subworkflows/local/drimseq_dexseq_dtu'
include { DRIMSEQ_DEXSEQ_DTU as DRIMSEQ_DEXSEQ_DTU_STAR_SALMON            } from '../subworkflows/local/drimseq_dexseq_dtu'
include { RMATS                                                           } from '../subworkflows/local/rmats'
include { DEXSEQ_DEU                                                      } from '../subworkflows/local/dexseq_deu'
include { EDGER_DEU                                                       } from '../subworkflows/local/edger_deu'
include { SUPPA as SUPPA_SALMON                                           } from '../subworkflows/local/suppa'
include { SUPPA as SUPPA_STAR_SALMON                                      } from '../subworkflows/local/suppa'
include { VISUALISE_MISO                                                  } from '../subworkflows/local/visualize_miso'
include { LEAFCUTTER                                                      } from '../subworkflows/local/leafcutter'

include { validateInputSamplesheet                                        } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { validateInputContrastsheet                                      } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { rmatsReadError                                                  } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { rmatsStrandednessError                                          } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { multiqcTsvFromList                                              } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'
include { paramsSummaryMultiqc                                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                          } from '../subworkflows/local/utils_nfcore_rnasplice_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList                                               } from 'plugin/nf-schema'
include { FASTQC                                                          } from '../modules/nf-core/fastqc/main'
include { SALMON_QUANT as SALMON_QUANT_SALMON                             } from '../modules/nf-core/salmon/quant'
include { SALMON_QUANT as SALMON_QUANT_STAR                               } from '../modules/nf-core/salmon/quant'
include { MULTIQC                                                         } from '../modules/nf-core/multiqc/'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                     } from '../modules/nf-core/custom/dumpsoftwareversions'
include { CAT_FASTQ                                                       } from '../modules/nf-core/cat/fastq'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE                                } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BAM_SORT_STATS_SAMTOOLS                                         } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { softwareVersionsToYAML                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                                                } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASPLICE {
    take:
    ch_samplesheet // channel: file(samplesheet)
    ch_contrastsheet  // channel: file(contrastsheet)
    ch_fasta // channel: path of genome fasta
    ch_gtf // channel: path of genome gtf
    ch_transcript_fasta // channel: path of transcript fasta
    ch_dexseq_gff
    ch_salmon_index // channel: path of salmon index
    ch_star_index // channel: path of star index
    ch_suppa_tpm
    ch_chrom_sizes
    is_aws_igenome
    multiqc_config // parameter value pass-through
    multiqc_logo // parameter value pass-through
    multiqc_methods_description // parameter value pass-through
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def pass_trimmed_reads = [:]

    //
    // Create channel from input file provided through params.input
    //
    if (params.source == "fastq") {
        ch_reads = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [meta.id, meta + [single_end: true], [fastq_1]]
                }
                else {
                    return [meta.id, meta + [single_end: false], [fastq_1, fastq_2]]
                }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map { meta, fastqs ->
                return [meta, fastqs.flatten()]
            }
    }
    else if (params.source == "genome_bam") {
        ch_reads = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input_genome_bam.json"))
            .map { meta, genome_bam ->
                def meta_map = [id: meta.id, condition: meta.condition]
                return [meta_map, [genome_bam]]
            }
    }
    else if (params.source == "transcriptome_bam") {
        ch_reads = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input_transcriptome_bam.json"))
            .map { meta, _genome_bam, transcriptome_bam ->
                def meta_map = [id: meta.id, condition: meta.condition]
                return [meta_map, [transcriptome_bam]]
            }
    }
    else if (params.source == "salmon_results") {
        ch_reads = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input_salmon_results.json"))
            .map { meta, salmon_results ->
                def meta_map = [id: meta.id, condition: meta.condition]
                return [meta_map, [salmon_results]]
            }
    }
    else {
        error("Invalid --source parameter: '${params.source}'. Must be one of: fastq, genome_bam, transcriptome_bam, salmon_results")
    }

    //
    // Create channel from contrasts file
    //
    channel.fromList(samplesheetToList(params.contrasts, "${projectDir}/assets/schema_contrasts.json"))
        .map { meta ->
            validateInputContrastsheet([[meta]])
            return [contrast: meta.contrast, treatment: meta.treatment, control: meta.control]
        }
        .collect()
        .flatMap { it }
        .set { ch_contrasts }

    // Branch samplesheet channel based on source type
    if (params.source == 'fastq') {
        ch_reads
            .map { meta, fastq ->
                def new_id = meta.id - ~/_T\d+/
                [meta + [id: new_id], fastq]
            }
            .groupTuple()
            .branch { meta, fastq ->
                single: fastq.size() == 1
                return [meta, fastq.flatten()]
                multiple: fastq.size() > 1
                return [meta, fastq.flatten()]
            }
            .set { ch_fastq }
    }
    else if (params.source == 'genome_bam') {
        ch_reads.set { ch_genome_bam }
    }
    else if (params.source == 'transcriptome_bam') {
        ch_reads.set { ch_transcriptome_bam }
    }
    else if (params.source == 'salmon_results') {
        ch_reads.set { ch_salmon_results }
    }

    // Check rMATS parameter configuration mapping checks
    if (params.rmats && params.source == 'fastq') {
        rmatsReadError(ch_reads)
        rmatsStrandednessError(ch_reads)
    }

    //
    // MODULE: Concatenate FastQ technical replicates if required
    //
    if (params.source == 'fastq') {
        CAT_FASTQ(ch_fastq.multiple)
        CAT_FASTQ.out.reads.mix(ch_fastq.single).set { ch_cat_fastq }
    }

    //
    // SUBWORKFLOW: Read QC and trim adapters with TrimGalore!
    //
    if (params.source == 'fastq') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
            ch_cat_fastq,
            params.skip_fastqc,
            false,
            true,
            params.skip_trimming,
            0,
            0,
        )
        ch_trim_reads = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_trim_read_count = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count

        // Parse custom evaluation warnings if reading maps fall below specific metric criteria
        ch_trim_read_count
            .map { meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return ["${meta.id}\t${num_reads}"]
                }
            }
            .collect()
            .map { tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_trimming_multiqc }
    }

    //
    // SUBWORKFLOW: BAM post-processing and indexing loops
    //
    if (params.source == 'genome_bam') {
        BAM_SORT_STATS_SAMTOOLS(ch_genome_bam, ch_fasta)
        ch_genome_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_genome_bam_index = BAM_SORT_STATS_SAMTOOLS.out.index
        ch_samtools_stats = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_samtools_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_samtools_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    }

    if (params.source == 'transcriptome_bam') {
        ch_transcriptome_bam_for_salmon = ch_transcriptome_bam
        BAM_SORT_STATS_SAMTOOLS(ch_transcriptome_bam, ch_fasta)
        ch_transcriptome_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_transcriptome_bam_index = BAM_SORT_STATS_SAMTOOLS.out.index
        ch_samtools_stats = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_samtools_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_samtools_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    }

    if (params.source == 'fastq' && !params.skip_alignment) {
        if (params.aligner == 'star' || params.aligner == 'star_salmon') {
            ALIGN_STAR(
                ch_trim_reads,
                ch_star_index.map { index -> [[:], index] },
                ch_gtf.map { gtf -> [[:], gtf] },
                params.star_ignore_sjdbgtf,
                '',
                params.seq_center ?: '',
                is_aws_igenome,
                ch_fasta.map { fasta -> [[:], fasta, []] },
            )
            ch_genome_bam = ALIGN_STAR.out.bam           // [meta, bam]
            ch_genome_bam_index = ALIGN_STAR.out.index
            ch_transcriptome_bam_for_salmon = ALIGN_STAR.out.bam_transcript
            ch_samtools_stats = ALIGN_STAR.out.stats
            ch_samtools_flagstat = ALIGN_STAR.out.flagstat
            ch_samtools_idxstats = ALIGN_STAR.out.idxstats
            ch_star_multiqc = ALIGN_STAR.out.log_final
        }
    }

    // ====================================================================
    // EXON SPLICING ENGINES
    // ====================================================================


    if (params.source == 'genome_bam' ||
            ( ( params.source == 'fastq') && ( !params.skip_alignment && (params.aligner == 'star' || params.aligner == 'star_salmon') ) )
        )
    {
        ch_dexseq_gff = params.gff_dexseq ?: ''

        if (params.dexseq_exon) {
            DEXSEQ_DEU(
                ch_gtf,
                ch_genome_bam,
                ch_dexseq_gff,
                ch_samplesheet,
                ch_contrastsheet,
                params.n_dexseq_plot,
                params.aggregation,
                params.alignment_quality,
            )
        }

        if (params.edger_exon) {
            EDGER_DEU(
                ch_gtf,
                ch_genome_bam,
                ch_samplesheet,
                ch_contrastsheet,
                params.n_edger_plot,
            )
        }

        if (params.rmats) {
            ch_genome_bam
                .map { meta, bam -> [meta.condition, meta, bam] }
                .set { ch_genome_bam_conditions }

            def lines
            if (params.input.startsWith('http')) {
                lines = new URL(params.input).text.readLines()
            }
            else {
                lines = new File(params.input).readLines()
            }
            def headers = lines[0].split(',')*.trim()
            def condition_idx = headers.indexOf('condition')
            def is_single_condition = lines[1..-1].collect { it -> it.split(',')[condition_idx].trim() }.unique().size() == 1

            RMATS(
                channel.value(file(params.input)),
                ch_contrastsheet,
                ch_genome_bam_conditions,
                ch_gtf,
                is_single_condition,
                params.rmats_read_len,
                params.rmats_splice_diff_cutoff,
                params.rmats_novel_splice_site,
                params.rmats_min_intron_len,
                params.rmats_max_exon_len,
                params.rmats_paired_stats,
            )
        }

        if (params.sashimi_plot == true) {
            VISUALISE_MISO(
                ch_gtf,
                ch_genome_bam,
                ch_genome_bam_index,
                params.fig_width,
                params.fig_height,
                params.miso_genes,
                params.miso_genes_file ?: false,
            )
        }

        if (params.leafcutter == true) {
            LEAFCUTTER(ch_genome_bam, ch_genome_bam_index)
        }
    }

    // ====================================================================
    // TRANSCRIPT ALIGNMENT & QUANTIFICATION ENGINES (SALMON)
    // ====================================================================

    if (params.source == 'transcriptome_bam' || ((params.source == 'fastq' && !params.skip_alignment && params.aligner == 'star_salmon'))) {
        ch_dummy_salmon_index = channel.value(file("${projectDir}/assets/dummy_file.txt", checkIfExists: true))

        SALMON_QUANT_STAR(
            ch_transcriptome_bam_for_salmon,
            ch_dummy_salmon_index,
            ch_gtf,
            ch_transcript_fasta,
            true,
            params.salmon_quant_libtype ?: '',
        )
        ch_salmon_results = SALMON_QUANT_STAR.out.results

        TX2GENE_TXIMPORT_STAR_SALMON(ch_salmon_results, ch_gtf)
        ch_tximport_tx2gene = TX2GENE_TXIMPORT_STAR_SALMON.out.tximport_tx2gene
        ch_txi_s = TX2GENE_TXIMPORT_STAR_SALMON.out.txi_s
        ch_txi_dtu = TX2GENE_TXIMPORT_STAR_SALMON.out.txi_dtu
        ch_txi_suppa_tpm = TX2GENE_TXIMPORT_STAR_SALMON.out.suppa_tpm

        if (params.dexseq_dtu) {
            ch_txi = (params.dtu_txi == "dtuScaledTPM") ? ch_txi_dtu : ch_txi_s
            DRIMSEQ_DEXSEQ_DTU_STAR_SALMON(
                ch_txi,
                ch_tximport_tx2gene,
                ch_samplesheet,
                ch_contrastsheet,
                params.min_samps_gene_expr,
                params.min_samps_feature_expr,
                params.min_samps_feature_prop,
                params.min_feature_expr,
                params.min_feature_prop,
                params.min_gene_expr,
            )
        }

        if (params.suppa) {
            ch_suppa_tpm = params.suppa_tpm ? ch_suppa_tpm : ch_txi_suppa_tpm
            SUPPA_STAR_SALMON(
                ch_gtf.map { gtf -> [[ id: gtf.baseName ], gtf] },
                ch_suppa_tpm.map { tpm_psi -> [[ id: tpm_psi.baseName ], tpm_psi] },
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
                params.diffsplice_isoform,
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
                params.clusterevents_isoform,
                params.clusterevents_dpsithreshold,
                params.clusterevents_eps,
                params.clusterevents_metric,
                params.clusterevents_min_pts,
                params.clusterevents_method,
                params.clusterevents_sigthreshold ?: false,
                params.clusterevents_separation ?: false,
                params.suppa_per_isoform,
            )
        }
    }

    if (params.source == 'fastq' && params.pseudo_aligner == 'salmon') {
        ch_dummy_transcript_fasta = channel.value(file("${projectDir}/assets/dummy_file2.txt", checkIfExists: true))

        SALMON_QUANT_SALMON(
            ch_trim_reads,
            ch_salmon_index,
            ch_gtf,
            ch_dummy_transcript_fasta,
            false,
            params.salmon_quant_libtype ?: '',
        )
        ch_salmon_results = SALMON_QUANT_SALMON.out.results

        TX2GENE_TXIMPORT_SALMON(ch_salmon_results, ch_gtf)
    }
    else if (params.source == 'salmon_results') {
        TX2GENE_TXIMPORT_SALMON(ch_salmon_results, ch_gtf)
    }

    if ((params.pseudo_aligner == 'salmon' && params.source == 'fastq') || (params.source == 'salmon_results')) {
        if (params.dexseq_dtu) {
            ch_txi = (params.dtu_txi == "dtuScaledTPM") ? TX2GENE_TXIMPORT_SALMON.out.txi_dtu : TX2GENE_TXIMPORT_SALMON.out.txi_s
            DRIMSEQ_DEXSEQ_DTU_SALMON(
                ch_txi,
                TX2GENE_TXIMPORT_SALMON.out.tximport_tx2gene,
                ch_samplesheet,
                ch_contrastsheet,
                params.min_samps_gene_expr,
                params.min_samps_feature_expr,
                params.min_samps_feature_prop,
                params.min_feature_expr,
                params.min_feature_prop,
                params.min_gene_expr,
            )
        }

        if (params.suppa) {
            ch_suppa_tpm = params.suppa_tpm ? ch_suppa_tpm : ch_txi_suppa_tpm
            SUPPA_SALMON(
                ch_gtf.map { gtf -> [ [ id: gtf.baseName ] , gtf] },
                ch_suppa_tpm.map { tpm_psi -> [ [ id: tpm_psi.baseName ] , tpm_psi] },
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
                params.diffsplice_isoform,
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
                params.clusterevents_isoform,
                params.clusterevents_dpsithreshold,
                params.clusterevents_eps,
                params.clusterevents_metric,
                params.clusterevents_min_pts,
                params.clusterevents_method,
                params.clusterevents_sigthreshold ?: false,
                params.clusterevents_separation ?: false,
                params.suppa_per_isoform,
            )
        }
    }

    if (params.isoformswitchanalyzer && (params.source == 'fastq' || params.source == 'salmon_results')) {
        ISOFORMSWITCHANALYZER(
            ch_salmon_results.collect { it -> it[1] },
            ch_gtf,
            ch_transcript_fasta,
            ch_samplesheet,
            ch_contrastsheet,
            params.isoformswitchanalyzer_alpha,
            params.isoformswitchanalyzer_dIF,
        )
    }

    //
    // Module: Generate stranded coverage bigWig files for visualisation in genome browsers
    //
    if ((params.source == 'genome_bam' && !params.skip_bigwig) || (!params.skip_alignment && !params.skip_bigwig && params.source == 'fastq')) {

        BEDTOOLS_GENOMECOV_FORWARD(
            ch_genome_bam.map { meta, bam -> [ meta, bam, 0 ]},
            ch_chrom_sizes.map { _meta, sizes -> sizes }, 'bedGraph', true)
        BEDTOOLS_GENOMECOV_REVERSE(ch_genome_bam.map { meta, bam -> [ meta, bam, 0]}, ch_chrom_sizes.map { _meta, sizes -> sizes }, 'bedGraph', true)

        BEDGRAPH_TO_BIGWIG_FORWARD(BEDTOOLS_GENOMECOV_FORWARD.out.genomecov, ch_chrom_sizes.map { _meta, sizes -> sizes })
        BEDGRAPH_TO_BIGWIG_REVERSE(BEDTOOLS_GENOMECOV_REVERSE.out.genomecov, ch_chrom_sizes.map { _meta, sizes -> sizes })
    }

    // ====================================================================
    // REPORTING LAYER & VERSION COLLATOR SYNC
    // ====================================================================

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnasplice_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    // Blend final channel maps
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    if (params.source == 'fastq') {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]))
    }

    if (params.pseudo_aligner == 'salmon' && params.source == 'fastq') {
        ch_multiqc_files = ch_multiqc_files.mix(ch_salmon_results.collect { it -> it[1] }.ifEmpty([]))
    }

    if (!params.skip_alignment && params.source == 'fastq' && (params.aligner == 'star_salmon' || params.aligner == "star")) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_multiqc.collect { it -> it[1] }.ifEmpty([]))
    }

    if ((params.source == 'genome_bam' || params.source == 'transcriptome_bam') || (!params.skip_alignment && params.source == 'fastq' && (params.aligner == 'star_salmon' || params.aligner == "star"))) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect { it -> it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect { it -> it[1] }.ifEmpty([]))

        if (params.edger_exon && params.source != 'transcriptome_bam') {
            ch_multiqc_files = ch_multiqc_files.mix(EDGER_DEU.out.featureCounts_summary.collect { it -> it[1] }.ifEmpty([]))
        }
    }

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'rnasplice'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
