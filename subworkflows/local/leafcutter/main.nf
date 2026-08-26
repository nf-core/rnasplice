include { REGTOOLS_JUNCTIONSEXTRACT      } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTERREGTOOLS     } from '../../../modules/nf-core/leafcutter/clusterregtools/main'

workflow LEAFCUTTER {

    take:
    ch_genome_bam
    ch_genome_bam_index

    main:

    REGTOOLS_JUNCTIONSEXTRACT(ch_genome_bam.join(ch_genome_bam_index), '')
    ch_juncs = REGTOOLS_JUNCTIONSEXTRACT.out.junc
        .map { _meta, junc -> junc }
        .collect()
        .map { junc_files -> [ [ id:'lc' ], junc_files ] }

    LEAFCUTTER_CLUSTERREGTOOLS(ch_juncs)

    emit:
    juncs                   = REGTOOLS_JUNCTIONSEXTRACT.out.junc.collect()
    counts                  = LEAFCUTTER_CLUSTERREGTOOLS.out.counts
    numers                  = LEAFCUTTER_CLUSTERREGTOOLS.out.numers
    pooled                  = LEAFCUTTER_CLUSTERREGTOOLS.out.pooled
    refined                 = LEAFCUTTER_CLUSTERREGTOOLS.out.refined
    sortedlibs              = LEAFCUTTER_CLUSTERREGTOOLS.out.sortedlibs
}
