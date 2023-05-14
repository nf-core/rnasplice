//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLIST as CREATE_BAMLIST_COND1  } from '../../modules/local/create_bamlist.nf'
include { CREATE_BAMLIST as CREATE_BAMLIST_COND2  } from '../../modules/local/create_bamlist.nf'
include { RMATS_PREP                              } from '../../modules/local/rmats_prep.nf'
include { RMATS_POST                              } from '../../modules/local/rmats_post.nf'
include { RMATS_PREP_SINGLE                       } from '../../modules/local/rmats_prep_single.nf'
include { RMATS_POST_SINGLE                       } from '../../modules/local/rmats_post_single.nf'

workflow RMATS {

    take:
    ch_contrastsheet          // channel: [ [contrast:A-B, treatment:A, control:B], [A1.meta, A2.meta], [B1.meta, B2.meta], [A1.bam A2.bam], [B1.bam B2.ba]]
    gtf                       // channel: /path/to/genome.gtf
    single_condition_bool     // channel: true/false

    main:

    ch_versions = Channel.empty()

    ch_contrastsheet_cond1 = ch_contrastsheet.map { [ it['treatment'], it['meta1'], it['bam1'] ] }

    ch_contrastsheet_cond2 = ch_contrastsheet.map { [ it['control'], it['meta2'], it['bam2'] ] }

    if (single_condition_bool) {

        //
        // Create input bam list file
        //

        CREATE_BAMLIST_COND1 (
            ch_contrastsheet_cond1
        )

        ch_bamlist_cond1 = CREATE_BAMLIST_COND1.out.bam_text.map{ it[1] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP_SINGLE (
            gtf,
            ch_bamlist_cond1,
            ch_contrastsheet_cond1
        )

        //
        // Run rMATs post step
        //

        RMATS_POST_SINGLE (
            gtf,
            ch_bamlist_cond1,
            ch_contrastsheet_cond1,
            RMATS_PREP_SINGLE.out.rmats_temp
        )

        ch_versions         = ch_versions.mix(RMATS_POST_SINGLE.out.versions)

        ch_rmats_prep       = RMATS_PREP_SINGLE.out.rmats_temp  //    path: rmats_temp/*
        ch_rmats_prep_log   = RMATS_PREP_SINGLE.out.log         //    path: rmats_prep.log
        ch_rmats_post       = RMATS_POST_SINGLE.out.rmats_post  //    path: rmats_post/*
        ch_rmats_post_log   = RMATS_POST_SINGLE.out.log         //    path: rmats_post.log

    } else {

        //
        // Create input bam list file
        //

        CREATE_BAMLIST_COND1 (
            ch_contrastsheet_cond1
        )

        CREATE_BAMLIST_COND2 (
            ch_contrastsheet_cond2
        )

        ch_bamlist_cond1 = CREATE_BAMLIST_COND1.out.bam_text.map{ it[1] }

        ch_bamlist_cond2 = CREATE_BAMLIST_COND2.out.bam_text.map{ it[1] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP (
            gtf,
            ch_bamlist_cond1,
            ch_bamlist_cond2,
            ch_contrastsheet_cond1,
            ch_contrastsheet_cond2
        )

        //
        // Run rMATs post step
        //

        RMATS_POST (
            gtf,
            ch_bamlist_cond1,
            ch_bamlist_cond2,
            ch_contrastsheet_cond1,
            ch_contrastsheet_cond2,
            RMATS_PREP.out.rmats_temp
        )

        ch_versions         = ch_versions.mix(RMATS_POST.out.versions)

        ch_rmats_prep       = RMATS_PREP.out.rmats_temp     //    path: rmats_prep/*
        ch_rmats_prep_log   = RMATS_PREP.out.log            //    path: rmats_prep.log
        ch_rmats_post       = RMATS_POST.out.rmats_post     //    path: rmats_post/*
        ch_rmats_post_log   = RMATS_POST.out.log            //    path: rmats_post.log

    }

    // Define output

    emit:

    rmats_prep       = ch_rmats_prep        //    path: rmats_temp/*
    rmats_prep_log   = ch_rmats_prep_log    //    path: rmats_prep.log
    rmats_post       = ch_rmats_post        //    path: rmats_post/*
    rmats_post_log   = ch_rmats_post_log    //    path: rmats_post.log

    versions = ch_versions.ifEmpty(null)    //    channel: [ versions.yml ]

}
