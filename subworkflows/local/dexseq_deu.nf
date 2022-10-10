//
// Dexseq DEU subworkflow
//

include { DEXSEQ_COUNT        } from '../../modules/local/dexseq_count'
include { DEXSEQ_EXON         } from '../../modules/local/dexseq_exon'

workflow DEXSEQ_DEU {

    take:

    ch_genome_bam_gff     // bam gff channel
    ch_dexseq_gff      // gff
    ch_samplesheet    // Channel.fromPath(params.input)
    read_method       // val: htseq or featurecounts       

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT (
        ch_genome_bam_gff
    )

    ch_versions = ch_versions.mix(DEXSEQ_COUNT.out.versions.first())

    //
    // MODULE: DEXSeq differential exon usage
    //

    DEXSEQ_EXON (
        DEXSEQ_COUNT.out.dexseq_clean_txt.collect(),
        ch_dexseq_gff,
        ch_samplesheet,
        read_method
    )

    ch_versions = ch_versions.mix(DEXSEQ_EXON.out.versions)

    emit:

    dexseq_clean_txt           = DEXSEQ_COUNT.out.dexseq_clean_txt.collect()  

    dexseq_exon_rds            = DEXSEQ_EXON.out.dexseq_exon_rds              // path: dxd.rds
    dexseq_exon_results_rds    = DEXSEQ_EXON.out.dexseq_exon_results_rds      // path: dxr.rds
    dexseq_exon_results_tsv    = DEXSEQ_EXON.out.dexseq_exon_results_tsv      // path: dxr.tsv
    qval_exon_rds              = DEXSEQ_EXON.out.qval_exon_rds                // path: qval.rds
    dexseq_exon_results_q_tsv  = DEXSEQ_EXON.out.dexseq_exon_results_q_tsv    // path: dxr.g.tsv

    versions                   = ch_versions                                  // channel: [ versions.yml ]
}