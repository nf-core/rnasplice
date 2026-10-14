//
// DEXSeq DEU subworkflow
//

include { DEXSEQ_PREPAREANNOTATION } from '../../../modules/nf-core/dexseq/prepareannotation'
include { DEXSEQ_COUNT } from '../../../modules/nf-core/dexseq/count'
include { DEXSEQ_EXON } from '../../../modules/nf-core/dexseq/exon'

workflow DEXSEQ_DEU {
    take:
    gtf // path gtf, only used when `skip_annotation` is false
    ch_genome_bam // bam channel
    ch_dexseq_gff // path dexseq gff, empty when `skip_annotation` is false
    skip_annotation // boolean: true when a DEXSeq annotation is supplied in `ch_dexseq_gff`
    ch_samplesheet // channel.fromPath(params.input)
    ch_contrastsheet // channel: [ contrast:, treatment:, control: ]

    main:

    if (!skip_annotation) {

        //
        // MODULE: DEXSeq Annotation
        //

        // `dexseq_prepare_annotation.py` iterates over Python sets, so a fixed hash seed
        // keeps the order of the aggregate gene IDs the same across runs
        DEXSEQ_PREPAREANNOTATION(
            gtf.map { gtf_file -> [[id: gtf_file.baseName], gtf_file] },
            0,
        )

        ch_dexseq_gff = DEXSEQ_PREPAREANNOTATION.out.gff
    }
    else {
        ch_dexseq_gff = ch_dexseq_gff.map { gff -> [[id: gff.baseName], gff] }
    }

    // The annotation is used by every sample and every contrast
    ch_dexseq_gff = ch_dexseq_gff.first()

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT(
        ch_genome_bam,
        ch_dexseq_gff,
    )

    //
    // MODULE: DEXSeq DEU
    //

    // Sort by file name so the collected count tables are emitted in a reproducible order.
    // `collect()` follows task completion order and `collect(sort: true)` sorts by the full
    // work directory path, neither of which is stable across runs. `DEXSEQ_EXON` itself does
    // not depend on this order: it builds the count file paths from the `sample` column of
    // the samplesheet
    ch_dexseq_clean_txt = DEXSEQ_COUNT.out.counts
        .map { _meta, txt -> txt }
        .toSortedList { txt_a, txt_b -> txt_a.name <=> txt_b.name }

    // `DEXSEQ_EXON` runs once per contrast, with the contrast name as its id
    ch_dexseq_exon_input = ch_contrastsheet
        .combine(ch_samplesheet)
        .combine(ch_dexseq_clean_txt.map { txts -> [txts] })
        .combine(ch_dexseq_gff)
        .multiMap { row, samplesheet, txts, gff_meta, gff ->
            def meta = [id: row.contrast]
            contrast: [meta, 'condition', row.control, row.treatment]
            counts: [meta, samplesheet, txts]
            gff: [gff_meta, gff]
        }

    DEXSEQ_EXON(
        ch_dexseq_exon_input.contrast,
        ch_dexseq_exon_input.counts,
        ch_dexseq_exon_input.gff,
    )

    emit:
    dexseq_clean_txt = ch_dexseq_clean_txt
    dexseq_exon_dataset_rds = DEXSEQ_EXON.out.dataset_rds
    dexseq_exon_results_rds = DEXSEQ_EXON.out.results_rds
    dexseq_gene_results_rds = DEXSEQ_EXON.out.gene_qvalues_rds
    dexseq_exon_results_tsv = DEXSEQ_EXON.out.results
    dexseq_gene_results_tsv = DEXSEQ_EXON.out.gene_qvalues
    dexseq_plot_results_pdf = DEXSEQ_EXON.out.plots
}
