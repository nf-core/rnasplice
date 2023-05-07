//
// Dexseq DTU subworkflow
//

include { DRIMSEQ_FILTER  } from '../../modules/local/drimseq_filter'
include { DEXSEQ_DTU      } from '../../modules/local/dexseq_dtu'
include { STAGER          } from '../../modules/local/stager'

workflow DRIMSEQ_DEXSEQ_DTU {

    take:

    txi                  // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    tximport_tx2gene     // path: tximport.tx2gene.tsv
    samplesheet          // path: /path/to/samplesheet.csv
    n_dexseq_plot        // val: numeric

    main:

    ch_versions = Channel.empty()

    //
    // Run DEXSeq DTU subworkflow - Filter, DEXSeq, and post-processing p-values with stageR.
    //

    DRIMSEQ_FILTER ( txi, tximport_tx2gene, samplesheet )

    ch_versions = ch_versions.mix(DRIMSEQ_FILTER.out.versions)

    DEXSEQ_DTU (
        DRIMSEQ_FILTER.out.drimseq_sample_data,
        DRIMSEQ_FILTER.out.drimseq_d_counts,
        n_dexseq_plot
    )

    ch_versions = ch_versions.mix(DEXSEQ_DTU.out.versions)

    def analysis_type = 'dexseq'

    DEXSEQ_DTU.out.dexseq_exon_results_rds.view()

    ch_dexseq_feature_rds = DEXSEQ_DTU.out.dexseq_exon_results_rds.collect { it ->
        [ it.baseName.toString().replaceAll("DEXSeqResults.", ""), it ]
    }

    ch_dexseq_gene_rds = DEXSEQ_DTU.out.dexseq_gene_results_rds.collect { it ->
        [ it.baseName.toString().replaceAll("perGeneQValue.", ""), it ]
    }

    ch_dexseq_feature_rds.view()
    ch_dexseq_gene_rds.view()

    ch_dexseq_results_rds = ch_dexseq_feature_rds.join(ch_dexseq_gene_rds)

    ch_dexseq_results_rds.view()

    STAGER (
        ch_dexseq_results_rds,
        analysis_type
    )

    ch_versions = ch_versions.mix(STAGER.out.versions)

    emit:

    drimseq_filter_rds    = DRIMSEQ_FILTER.out.drimseq_filter_rds
    drimseq_sample_data   = DRIMSEQ_FILTER.out.drimseq_sample_data
    drimseq_d_counts      = DRIMSEQ_FILTER.out.drimseq_d_counts

    dexseq_exon_dataset_rds  = DEXSEQ_DTU.out.dexseq_exon_dataset_rds
    dexseq_exon_results_rds  = DEXSEQ_DTU.out.dexseq_exon_results_rds
    dexseq_exon_results_tsv  = DEXSEQ_DTU.out.dexseq_exon_results_tsv
    dexseq_gene_results_rds  = DEXSEQ_DTU.out.dexseq_gene_results_rds
    dexseq_gene_results_tsv  = DEXSEQ_DTU.out.dexseq_gene_results_tsv

    stager_rds            = STAGER.out.stager_rds
    stager_padj_rds       = STAGER.out.stager_padj_rds
    stager_padj_tsv       = STAGER.out.stager_padj_tsv

    versions              = ch_versions                                  // channel: [ versions.yml ]
}
