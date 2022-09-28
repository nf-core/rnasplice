//
// Dexseq DTU subworkflow
//

include { DRIMSEQ_FILTER  } from '../../modules/local/drimseq_filter'
include { DEXSEQ_DTU         } from '../../modules/local/dexseq_dtu'
include { STAGER         } from '../../modules/local/stager'

workflow DRIMSEQ_DEXSEQ_DTU {

    take:

    txi                  // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    tximport_tx2gene     // path: tximport.tx2gene.tsv
    samplesheet          // path: /path/to/samplesheet.csv

    main:

    ch_versions = Channel.empty()

    //
    // Run DEXSeq DTU subworkflow - Filter, DEXSeq, and post-processing p-values with stageR.
    //

    DRIMSEQ_FILTER ( txi, tximport_tx2gene, samplesheet )

    ch_versions = ch_versions.mix(DRIMSEQ_FILTER.out.versions)

    DEXSEQ_DTU ( DRIMSEQ_FILTER.out.drimseq_filter_rds )

    ch_versions = ch_versions.mix(DEXSEQ_DTU.out.versions)

    def analysis_type = 'dexseq'
    
    STAGER ( 
        DEXSEQ_DTU.out.dexseq_results_rds,
        analysis_type,
        DEXSEQ_DTU.out.qval_rds
    )

    ch_versions = ch_versions.mix(STAGER.out.versions)

    emit:

    drimseq_filter_rds    = DRIMSEQ_FILTER.out.drimseq_filter_rds        // path: d.rds

    dexseq_rds            = DEXSEQ_DTU.out.dexseq_rds                    // path: dxd.rds
    dexseq_results_rds    = DEXSEQ_DTU.out.dexseq_results_rds            // path: dxr.rds
    dexseq_results_tsv    = DEXSEQ_DTU.out.dexseq_results_tsv            // path: dxr.tsv
    qval_rds              = DEXSEQ_DTU.out.qval_rds                      // path: qval.rds
    dexseq_results_q_tsv  = DEXSEQ_DTU.out.dexseq_results_q_tsv          // path: dxr.g.tsv

    stager_padj_tsv       = STAGER.out.stager_padj_tsv                   // path: *.stageR.padj.tsv
    stager_padj_rds       = STAGER.out.stager_padj_rds                   // path: **.stageR.padj.rds
    stager_rds            = STAGER.out.stager_rds                        // path: *.stageRObj.rds

    versions              = ch_versions                                  // channel: [ versions.yml ]
}