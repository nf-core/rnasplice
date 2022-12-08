//
// edgeR DEU subworkflow
//

include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_EDGER } from '../../modules/nf-core/subread/featurecounts/main'
include { EDGER_EXON            } from '../../modules/local/edger_exon'

workflow EDGER_DEU {

    take:

    gtf                  // path: gtf
    ch_genome_bam        // channel: [ val(meta), path(bams) ]
    ch_samplesheet       // Channel.fromPath(params.input)

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: EDGER_SUBREAD_FEATURECOUNTS
    //

    ch_feature_counts = ch_genome_bam.combine(gtf)

    SUBREAD_FEATURECOUNTS_EDGER(
        ch_feature_counts
    )

    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_EDGER.out.versions.first())

    //
    // MODULE: EDGER_EXON
    //
    EDGER_EXON (
        SUBREAD_FEATURECOUNTS_EDGER.out.counts.collect({it[1]}),
        ch_samplesheet
    )

    emit:

    featureCounts_summary  = SUBREAD_FEATURECOUNTS_EDGER.out.summary            // path featureCounts.txt.summary

    versions               = ch_versions                        // channel: [ versions.yml ]

}
