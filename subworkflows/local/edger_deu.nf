//
// edgeR DEU subworkflow
//

include { SUBREAD_FLATTENGTF    } from '../../modules/local/flattengtf'
include { SUBREAD_FEATURECOUNTS } from '../../modules/nf-core/subread/featurecounts/main'
include { EDGER_EXON            } from '../../modules/local/edger_exon'

workflow EDGER_DEU {

    take:

    gtf                  // path: gtf
    ch_genome_bam        // channel: [ val(meta), path(bams) ]
    ch_samplesheet       // Channel.fromPath(params.input)
    ch_contrastsheet     // Channel.fromPath(params.contrasts)
    n_edger_plot         // val: integer to plot

    main:

    ch_versions = Channel.empty()


    // MODULE: SUBREAD_FLATTENGTF

    SUBREAD_FLATTENGTF(gtf)

    //
    // MODULE: SUBREAD_FEATURECOUNTS
    //

    ch_feature_counts = ch_genome_bam.combine(SUBREAD_FLATTENGTF.out.saf)

    SUBREAD_FEATURECOUNTS(ch_feature_counts)

    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // MODULE: EDGER_COUNTS AND PLOT
    //
    EDGER_EXON (
        SUBREAD_FEATURECOUNTS.out.counts.collect({it[1]}),
        ch_samplesheet,
        ch_contrastsheet,
        n_edger_plot
    )


    emit:

    featureCounts_summary  = SUBREAD_FEATURECOUNTS.out.summary  // path featureCounts.txt.summary
    versions               = ch_versions                        // channel: [ versions.yml ]

}
