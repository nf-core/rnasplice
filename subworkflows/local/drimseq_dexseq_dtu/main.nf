//
// DRIMSeq filtering and DEXSeq DTU subworkflow
//

include { DRIMSEQ_DMFILTER } from '../../../modules/local/drimseq/dmfilter'
include { DEXSEQ_DTU } from '../../../modules/nf-core/dexseq/dtu'
include { STAGER } from '../../../modules/local/stager'
include { GENEANNOTATION_ANNOTATE as GENEANNOTATION_ANNOTATE_DEXSEQ } from '../../../modules/local/geneannotation/annotate'
include { GENEANNOTATION_ANNOTATE as GENEANNOTATION_ANNOTATE_STAGER } from '../../../modules/local/geneannotation/annotate'

workflow DRIMSEQ_DEXSEQ_DTU {
    take:
    txi // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    tximport_tx2gene // path: tximport.tx2gene.tsv
    samplesheet // path: /path/to/samplesheet.csv
    ch_contrastsheet // channel: [ contrast:, treatment:, control:, blocking: (when set) ]
    ch_gene_annotation // channel: [ val(meta), path(tsv) ], the gene annotation table, as a value channel
    stager_alpha // float: target overall false discovery rate of stageR

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

    // `DEXSEQ_DTU` runs once per contrast, with the contrast name as its id and, when
    // set, its blocking variables, which `ext.args` passes to the module. DRIMSeq
    // filters features, not samples, so the samplesheet still describes every column of
    // the filtered count table
    ch_dexseq_dtu_input = ch_contrastsheet
        .combine(samplesheet)
        .combine(DRIMSEQ_DMFILTER.out.drimseq_counts_tsv)
        .multiMap { row, samplesheet_csv, counts_tsv ->
            def meta = row.blocking ? [id: row.contrast, blocking: row.blocking] : [id: row.contrast]
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
        stager_alpha,
    )

    //
    // MODULE: Copies of the DEXSeq and stageR tables with the gene names
    //

    GENEANNOTATION_ANNOTATE_DEXSEQ(
        DEXSEQ_DTU.out.results.mix(DEXSEQ_DTU.out.gene_qvalues),
        ch_gene_annotation,
        'groupID',
    )

    GENEANNOTATION_ANNOTATE_STAGER(
        STAGER.out.stager_padj_tsv,
        ch_gene_annotation,
        'geneID',
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
    dexseq_annotated_tsv = GENEANNOTATION_ANNOTATE_DEXSEQ.out.annotated // channel: [ val(meta), path(tsv) ], the transcript results and gene q-values with the gene names
    stager_padj_annotated_tsv = GENEANNOTATION_ANNOTATE_STAGER.out.annotated // channel: [ val(meta), path(tsv) ]
}
