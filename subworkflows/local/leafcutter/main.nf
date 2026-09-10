//
// LeafCutter intron clustering subworkflow
//

include { REGTOOLS_JUNCTIONSEXTRACT } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTERREGTOOLS } from '../../../modules/nf-core/leafcutter/clusterregtools/main'

workflow LEAFCUTTER {
    take:
    ch_genome_bam // channel: [ val(meta), path(bam) ]
    ch_genome_bam_index // channel: [ val(meta), path(bai) ]

    main:

    //
    // MODULE: REGTOOLS_JUNCTIONSEXTRACT
    //

    REGTOOLS_JUNCTIONSEXTRACT(ch_genome_bam.join(ch_genome_bam_index), '')

    //
    // MODULE: LEAFCUTTER_CLUSTERREGTOOLS
    //

    // The order the junction files are collected in does not reach the results:
    // `LEAFCUTTER_CLUSTERREGTOOLS` sorts the file list itself before handing it to
    // `leafcutter-cluster`, so no sorting is needed here
    ch_juncs = REGTOOLS_JUNCTIONSEXTRACT.out.junc
        .map { _meta, junc -> junc }
        .collect()
        .map { junc_files -> [[id: 'lc'], junc_files] }

    LEAFCUTTER_CLUSTERREGTOOLS(ch_juncs)

    emit:
    juncs = REGTOOLS_JUNCTIONSEXTRACT.out.junc // channel: [ val(meta), path(junc) ]
    counts = LEAFCUTTER_CLUSTERREGTOOLS.out.counts // channel: [ val(meta), path(counts) ]
    numers = LEAFCUTTER_CLUSTERREGTOOLS.out.numers // channel: [ val(meta), path(numers) ]
    pooled = LEAFCUTTER_CLUSTERREGTOOLS.out.pooled // channel: [ val(meta), path(pooled) ]
    refined = LEAFCUTTER_CLUSTERREGTOOLS.out.refined // channel: [ val(meta), path(refined) ]
    sortedlibs = LEAFCUTTER_CLUSTERREGTOOLS.out.sortedlibs // channel: [ val(meta), path(sortedlibs) ]
}
