//
// DRIMSeq filtering and DEXSeq DTU subworkflow
//

include { DRIMSEQ_DMFILTER } from '../../../modules/local/drimseq/dmfilter'
include { DEXSEQ_DTU } from '../../../modules/nf-core/dexseq/dtu'
include { STAGER } from '../../../modules/local/stager'

workflow DRIMSEQ_DEXSEQ_DTU {
    take:
    txi // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    tximport_tx2gene // path: tximport.tx2gene.tsv
    samplesheet // path: /path/to/samplesheet.csv
    ch_contrastsheet // channel: [ contrast:, treatment:, control: ]

    main:

    //
    // MODULE: DRIMSeq filter
    //

    DRIMSEQ_DMFILTER(
        txi,
        tximport_tx2gene,
        samplesheet,
    )

    //
    // MODULE: DEXSeq DTU
    //

    // `DEXSEQ_DTU` runs once per contrast, with the contrast name as its id. DRIMSeq
    // filters features, not samples, so the samplesheet still describes every column of
    // the filtered count table
    ch_dexseq_dtu_input = ch_contrastsheet
        .combine(samplesheet)
        .combine(DRIMSEQ_DMFILTER.out.drimseq_counts_tsv)
        .multiMap { row, samplesheet_csv, counts_tsv ->
            def meta = [id: row.contrast]
            contrast: [meta, 'condition', row.control, row.treatment]
            counts: [meta, samplesheet_csv, counts_tsv]
        }

    DEXSEQ_DTU(
        ch_dexseq_dtu_input.contrast,
        ch_dexseq_dtu_input.counts,
    )

    //
    // MODULE: stageR
    //

    ch_dexseq_feature_gene_tsv = DEXSEQ_DTU.out.results.join(DEXSEQ_DTU.out.gene_qvalues)

    def analysis_type = 'dexseq'

    STAGER(
        ch_dexseq_feature_gene_tsv,
        analysis_type,
    )

    emit:
    drimseq_dataset_rds = DRIMSEQ_DMFILTER.out.drimseq_dataset_rds
    drimseq_samples_tsv = DRIMSEQ_DMFILTER.out.drimseq_samples_tsv
    drimseq_counts_tsv = DRIMSEQ_DMFILTER.out.drimseq_counts_tsv
    dexseq_exon_dataset_rds = DEXSEQ_DTU.out.dataset_rds
    dexseq_exon_results_rds = DEXSEQ_DTU.out.results_rds
    dexseq_exon_results_tsv = DEXSEQ_DTU.out.results
    dexseq_gene_results_rds = DEXSEQ_DTU.out.gene_qvalues_rds
    dexseq_gene_results_tsv = DEXSEQ_DTU.out.gene_qvalues
    stager_rds = STAGER.out.stager_rds
    stager_padj_rds = STAGER.out.stager_padj_rds
    stager_padj_tsv = STAGER.out.stager_padj_tsv
}
