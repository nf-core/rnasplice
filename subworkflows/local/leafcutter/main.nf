//
// LeafCutter intron clustering subworkflow
//

include { REGTOOLS_JUNCTIONSEXTRACT as REGTOOLS_JUNCTIONSEXTRACT_FR } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { REGTOOLS_JUNCTIONSEXTRACT as REGTOOLS_JUNCTIONSEXTRACT_RF } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { REGTOOLS_JUNCTIONSEXTRACT as REGTOOLS_JUNCTIONSEXTRACT_XS } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTERREGTOOLS } from '../../../modules/nf-core/leafcutter/clusterregtools/main'

workflow LEAFCUTTER {
    take:
    ch_genome_bam // channel: [ val(meta), path(bam) ]
    ch_genome_bam_index // channel: [ val(meta), path(bai) ]

    main:

    //
    // MODULE: REGTOOLS_JUNCTIONSEXTRACT
    //

    // `leafcutter-cluster` discards every junction that has no strand, so regtools has to
    // be told where to take the strand from. A stranded library has it in the orientation
    // of its reads, which regtools reads off the BAM: `FR` for a forward (second strand)
    // library, `RF` for a reverse (first strand) one. An unstranded library has no
    // orientation to read it from and is left with the aligner `XS` tag, which STAR only
    // writes when asked to infer the strand from the splice motif, see `conf/modules.config`
    ch_bam_bai = ch_genome_bam
        .join(ch_genome_bam_index)
        .branch { meta, _bam, _bai ->
            forward: meta.strandedness == 'forward'
            reverse: meta.strandedness == 'reverse'
            unstranded: true
        }

    REGTOOLS_JUNCTIONSEXTRACT_FR(ch_bam_bai.forward, 'FR')
    REGTOOLS_JUNCTIONSEXTRACT_RF(ch_bam_bai.reverse, 'RF')
    REGTOOLS_JUNCTIONSEXTRACT_XS(ch_bam_bai.unstranded, 'XS')

    ch_junc = REGTOOLS_JUNCTIONSEXTRACT_FR.out.junc
        .mix(REGTOOLS_JUNCTIONSEXTRACT_RF.out.junc)
        .mix(REGTOOLS_JUNCTIONSEXTRACT_XS.out.junc)

    //
    // MODULE: LEAFCUTTER_CLUSTERREGTOOLS
    //

    // The order the junction files are collected in does not reach the results:
    // `LEAFCUTTER_CLUSTERREGTOOLS` sorts the file list itself before handing it to
    // `leafcutter-cluster`, so no sorting is needed here
    ch_juncs = ch_junc
        .map { _meta, junc -> junc }
        .collect()
        .map { junc_files -> [[id: 'lc'], junc_files] }

    LEAFCUTTER_CLUSTERREGTOOLS(ch_juncs)

    emit:
    juncs = ch_junc // channel: [ val(meta), path(junc) ]
    counts = LEAFCUTTER_CLUSTERREGTOOLS.out.counts // channel: [ val(meta), path(counts) ]
    numers = LEAFCUTTER_CLUSTERREGTOOLS.out.numers // channel: [ val(meta), path(numers) ]
    pooled = LEAFCUTTER_CLUSTERREGTOOLS.out.pooled // channel: [ val(meta), path(pooled) ]
    refined = LEAFCUTTER_CLUSTERREGTOOLS.out.refined // channel: [ val(meta), path(refined) ]
    sortedlibs = LEAFCUTTER_CLUSTERREGTOOLS.out.sortedlibs // channel: [ val(meta), path(sortedlibs) ]
}
