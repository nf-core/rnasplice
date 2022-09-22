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

    STAGER ( DEXSEQ.out.dexseq_rds )

    ch_versions = ch_versions.mix(STAGER.out.versions)

    emit:

    drimseq_filter_csv    = DRIMSEQ_FILTER.out.drimseq_filter_csv    // path: *.csv
    drimseq_filter_rds    = DRIMSEQ_FILTER.out.drimseq_filter_rds    // path: *.rds

    dexseq_csv            = DEXSEQ.out.dexseq_csv                    // path: *.csv
    dexseq_rds            = DEXSEQ.out.dexseq_rds                    // path: *.rds

    stager_csv            = STAGER.out.stager_csv                    // path: *.csv
    stager_rds            = STAGER.out.stager_rds                    // path: *.rds

    versions              = ch_versions                              // channel: [ versions.yml ]
}