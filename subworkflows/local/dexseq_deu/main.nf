//
// DEXSeq DEU subworkflow
//

include { DEXSEQ_ANNOTATION } from '../../../modules/local/dexseq/annotation'
include { DEXSEQ_COUNT } from '../../../modules/local/dexseq/count'
include { DEXSEQ_EXON } from '../../../modules/local/dexseq/exon'

workflow DEXSEQ_DEU {
    take:
    gtf // path gtf
    ch_genome_bam // bam channel
    ch_dexseq_gff // path dexseq gff
    ch_samplesheet // channel.fromPath(params.input)
    ch_contrastsheet // channel.fromPath()
    n_dexseq_plot // val: numeric
    aggregation // params.aggregation
    alignment_quality // params.alignment_quality

    main:

    if (!params.gff_dexseq) {

        //
        // MODULE: DEXSeq Annotation
        //

        DEXSEQ_ANNOTATION(
            gtf,
            aggregation,
        )

        ch_dexseq_gff = DEXSEQ_ANNOTATION.out.gff
    }

    ch_genome_bam_gff = ch_genome_bam.combine(ch_dexseq_gff)

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT(
        ch_genome_bam_gff,
        alignment_quality,
    )

    //
    // MODULE: DEXSeq DEU
    //

    // Sort by file name so the collected count tables are emitted in a reproducible order.
    // `collect()` follows task completion order and `collect(sort: true)` sorts by the full
    // work directory path, neither of which is stable across runs. `DEXSEQ_EXON` itself does
    // not depend on this order: `run_dexseq_exon.R` builds the count file paths from the
    // `sample` column of the samplesheet
    ch_dexseq_clean_txt = DEXSEQ_COUNT.out.dexseq_clean_txt
        .map { _meta, txt -> txt }
        .toSortedList { txt_a, txt_b -> txt_a.name <=> txt_b.name }

    DEXSEQ_EXON(
        ch_dexseq_clean_txt,
        ch_dexseq_gff,
        ch_samplesheet,
        ch_contrastsheet,
        n_dexseq_plot,
    )

    emit:
    dexseq_clean_txt = ch_dexseq_clean_txt
    dexseq_exon_dataset_rds = DEXSEQ_EXON.out.dexseq_exon_dataset_rds
    dexseq_exon_results_rds = DEXSEQ_EXON.out.dexseq_exon_results_rds
    dexseq_gene_results_rds = DEXSEQ_EXON.out.dexseq_gene_results_rds
    dexseq_exon_results_csv = DEXSEQ_EXON.out.dexseq_exon_results_csv
    dexseq_gene_results_csv = DEXSEQ_EXON.out.dexseq_gene_results_csv
    dexseq_plot_results_pdf = DEXSEQ_EXON.out.dexseq_plot_results_pdf
}
