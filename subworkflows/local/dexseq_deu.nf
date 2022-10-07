//
// Dexseq DEU subworkflow
//

include { DEXSEQ_ANNOTATION   } from '../../modules/local/dexseq_annotation'
include { DEXSEQ_COUNT        } from '../../modules/local/dexseq_count'
include { DEXSEQ_EXON         } from '../../modules/local/dexseq_exon'

workflow DEXSEQ_DEU {

    take:

    gtf               // path: gtf 
    ch_genome_bam     // bam channel
    ch_samplesheet    // Channel.fromPath(params.input)
    read_method       // val: htseq or featurecounts

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Preparing the annotation using DEXSeq
    //

    DEXSEQ_ANNOTATION (
        gtf
    )

    ch_versions = ch_versions.mix(DEXSEQ_ANNOTATION.out.versions.first())

    //
    // MODULE: DEXSeq Count
    //

    DEXSEQ_COUNT (
        DEXSEQ_ANNOTATION.out.gff,
        ch_genome_bam
    )

    ch_versions = ch_versions.mix(DEXSEQ_COUNT.out.versions.first())

    //
    // MODULE: DEXSeq differential exon usage
    //

    // ch_samplesheet = Channel.fromPath(params.input)
    // def read_method = "htseq"

    DEXSEQ_EXON (
        DEXSEQ_COUNT.out.dexseq_clean_txt.collect(),
        DEXSEQ_ANNOTATION.out.gff,
        ch_samplesheet,
        read_method
    )

    ch_versions = ch_versions.mix(DEXSEQ_EXON.out.versions)

    emit:

    gff                        = DEXSEQ_ANNOTATION.out.gff                    // path: gff

    dexseq_clean_txt           = DEXSEQ_COUNT.out.dexseq_clean_txt.collect()  

    dexseq_exon_rds            = DEXSEQ_EXON.out.dexseq_exon_rds              // path: dxd.rds
    dexseq_exon_results_rds    = DEXSEQ_EXON.out.dexseq_exon_results_rds      // path: dxr.rds
    dexseq_exon_results_tsv    = DEXSEQ_EXON.out.dexseq_exon_results_tsv      // path: dxr.tsv
    qval_exon_rds              = DEXSEQ_EXON.out.qval_exon_rds                // path: qval.rds
    dexseq_exon_results_q_tsv  = DEXSEQ_EXON.out.dexseq_exon_results_q_tsv    // path: dxr.g.tsv

    versions                   = ch_versions                                  // channel: [ versions.yml ]
}