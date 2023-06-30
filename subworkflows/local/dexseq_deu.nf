//
// DEXSeq DEU subworkflow
//

include { DEXSEQ_ANNOTATION   } from '../../modules/local/dexseq_annotation'
include { DEXSEQ_COUNT        } from '../../modules/local/dexseq_count'
include { DEXSEQ_EXON         } from '../../modules/local/dexseq_exon'

workflow DEXSEQ_DEU {

    take:

    gtf                // path gtf
    ch_genome_bam      // bam channel
    ch_dexseq_gff      // path dexseq gff
    ch_samplesheet     // Channel.fromPath(params.input)
    ch_contrastsheet   // Channel.fromPath()
    n_dexseq_plot      // val: numeric
    aggregation        // params.aggregation
    alignment_quality  // params.alignment_quality

    main:

    ch_versions = Channel.empty()

    if (!ch_dexseq_gff) {

        //
        // MODULE: DEXSeq Annotation
        //

        DEXSEQ_ANNOTATION (
            gtf,
            aggregation
        )

        ch_versions = ch_versions.mix(DEXSEQ_ANNOTATION.out.versions.first())

        ch_dexseq_gff = DEXSEQ_ANNOTATION.out.gff

    }

    ch_genome_bam_gff = ch_genome_bam.combine(ch_dexseq_gff)

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT (
        ch_genome_bam_gff,
        alignment_quality
    )

    ch_versions = ch_versions.mix(DEXSEQ_COUNT.out.versions.first())

    //
    // MODULE: DEXSeq DEU
    //

    DEXSEQ_EXON (
        DEXSEQ_COUNT.out.dexseq_clean_txt.map{ it[1] }.collect(),
        ch_dexseq_gff,
        ch_samplesheet,
        ch_contrastsheet,
        n_dexseq_plot
    )

    ch_versions = ch_versions.mix(DEXSEQ_EXON.out.versions)

    emit:

    dexseq_clean_txt        = DEXSEQ_COUNT.out.dexseq_clean_txt.map{ it[1] }.collect()

    dexseq_exon_dataset_rds = DEXSEQ_EXON.out.dexseq_exon_dataset_rds
    dexseq_exon_results_rds = DEXSEQ_EXON.out.dexseq_exon_results_rds
    dexseq_gene_results_rds = DEXSEQ_EXON.out.dexseq_gene_results_rds
    dexseq_exon_results_csv = DEXSEQ_EXON.out.dexseq_exon_results_csv
    dexseq_gene_results_csv = DEXSEQ_EXON.out.dexseq_gene_results_csv
    dexseq_plot_results_pdf = DEXSEQ_EXON.out.dexseq_plot_results_pdf

    versions                = ch_versions                                  // channel: [ versions.yml ]

}
