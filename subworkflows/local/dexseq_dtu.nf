//
// Dexseq DTU subworkflow
//

include { DRIMSEQ_FILTER  } from '../../modules/local/drimseq_filter'
include { DEXSEQ         } from '../../modules/local/dexseq'
include { STAGER         } from '../../modules/local/stager'

workflow DEXSEQ_DTU {

    take:

    txi         // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    tx2gene     // path: *.tx2gene.tsv
    samplesheet // path: /path/to/samplesheet.csv

    main:

    ch_versions = Channel.empty()

    //
    // Run DEXSeq DTU subworkflow - Filter, DEXSeq, and post-processing p-values with stageR.
    //

    DRIMSEQ_FILTER ( txi, tx2gene, samplesheet )

    ch_versions = ch_versions.mix(DRIMSEQ_FILTER.out.versions)

    DEXSEQ ( DRIMSEQ_FILTER.out.drimseq_filter_rds )

    ch_versions = ch_versions.mix(DEXSEQ.out.versions)

    def analysis_type = 'dexseq'
    
    STAGER ( 
        DEXSEQ.out.dexseq_rds,
        analysis_type,
        tx2gene,
        DEXSEQ.out.qval_rds
    )

    ch_versions = ch_versions.mix(STAGER.out.versions)

    emit:

    drimseq_filter_rds    = DRIMSEQ_FILTER.out.drimseq_filter_rds    // path: d.rds

    dexseq_rds            = DEXSEQ.out.dexseq_rds                    // path: dxd.rds
    dexseq_results_rds    = DEXSEQ.out.dexseq_results_rds            // path: dxr.rds
    dexseq_results_tsv    = DEXSEQ.out.dexseq_results_tsv            // path: dxr.tsv
    qval_rds              = DEXSEQ.out.qval_rds                      // path: qval.rds
    dexseq_results_q_tsv  = DEXSEQ.out.dexseq_results_q_tsv          // path: dxr.g.tsv

    stager_padj_tsv       = STAGER.out.stager_padj_tsv               // path: *.stageR.padj.tsv
    stager_padj_rds       = STAGER.out.stager_padj_rds               // path: **.stageR.padj.rds
    stager_rds            = STAGER.out.stager_rds                    // path: *.stageRObj.rds

    versions              = ch_versions                              // channel: [ versions.yml ]
}