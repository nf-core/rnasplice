//
// edgeR DEU subworkflow
//

include { SUBREAD_FEATURECOUNTS } from '../../modules/nf-core/subread/featurecounts/main'
include { EDGER_EXON            } from '../../modules/local/edger_exon'

workflow EDGER_DEU {

    take:
    
    gtf                  // path: gtf
    ch_genome_bam        // channel: [ val(meta), path(bams) ]
    ch_samplesheet       // Channel.fromPath(params.input)

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: SUBREAD_FEATURECOUNTS
    //

    ch_feature_counts = ch_genome_bam.combine(gtf)

    SUBREAD_FEATURECOUNTS(ch_feature_counts)

    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // MODULE: EDGER_COUNTS
    //
    EDGER_EXON (
        SUBREAD_FEATURECOUNTS.out.counts.collect({it[1]}),
        ch_samplesheet
    )

    //
    // MODULE: EDGER_EXON
    //

    emit:
    
    versions = ch_versions      // channel: [ versions.yml ]

}