//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLIST         } from '../../modules/local/create_bamlist.nf'
include { RMATS_PREP             } from '../../modules/local/rmats_prep.nf'
include { RMATS_POST             } from '../../modules/local/rmats_post.nf'

include { CREATE_BAMLIST_SINGLE  } from '../../modules/local/create_bamlist_single.nf'
include { RMATS_PREP_SINGLE      } from '../../modules/local/rmats_prep_single.nf'
include { RMATS_POST_SINGLE      } from '../../modules/local/rmats_post_single.nf'

workflow RMATS {

    take:
    ch_contrast            // channel: [ contrast:A-B, treatment:A, control:B, meta1:[A1.meta, A2.meta], meta2:[B1.meta, B2.meta], bam1:[A1.bam A2.bam], bam2:[B1.bam B2.ba]]
    gtf                    // channel: /path/to/genome.gtf
    single_condition_bool  // channel: true/false

    main:

    ch_versions = Channel.empty()

    if (single_condition_bool) {

        //
        // Create input channels
        //

        ch_contrast_cond1 = ch_contrast.map { [ it.treatment, it.meta1, it.bam1 ] }

        //
        // Create input bam list file
        //

        CREATE_BAMLIST_SINGLE (
            ch_contrast_cond1
        )

        //
        // Create input channels
        //

        ch_bamlist_cond1 = CREATE_BAMLIST_SINGLE.out.bam_text.map { it[3] }

        ch_contrast_cond1 = CREATE_BAMLIST_SINGLE.out.bam_text.map { [ it[0], it[1], it[2] ] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP_SINGLE (
            gtf,
            ch_bamlist_cond1,
            ch_contrast_cond1
        )

        //
        // Run rMATs post step
        //

        RMATS_POST_SINGLE (
            gtf,
            ch_bamlist_cond1,
            ch_contrast_cond1,
            RMATS_PREP_SINGLE.out.rmats_temp
        )

        ch_versions         = ch_versions.mix(RMATS_POST_SINGLE.out.versions)

        ch_rmats_prep       = RMATS_PREP_SINGLE.out.rmats_temp  //    path: rmats_temp/*
        ch_rmats_prep_log   = RMATS_PREP_SINGLE.out.log         //    path: rmats_prep.log
        ch_rmats_post       = RMATS_POST_SINGLE.out.rmats_post  //    path: rmats_post/*
        ch_rmats_post_log   = RMATS_POST_SINGLE.out.log         //    path: rmats_post.log

    } else {

        //
        // Create input channels
        //

        ch_contrast_cond1 = ch_contrast.map { [ it.treatment, it.meta1, it.bam1 ] }

        ch_contrast_cond2 = ch_contrast.map { [ it.control, it.meta2, it.bam2 ] }

        //
        // Create input bam list file
        //

        CREATE_BAMLIST (
            ch_contrast_cond1,
            ch_contrast_cond2
        )

        //
        // Create input channels
        //

        ch_bamlist_cond1 = CREATE_BAMLIST.out.bam1_text.map { it[3] }

        ch_contrast_cond1 = CREATE_BAMLIST.out.bam1_text.map { [ it[0], it[1], it[2] ] }

        ch_bamlist_cond2 = CREATE_BAMLIST.out.bam2_text.map { it[3] }

        ch_contrast_cond2 = CREATE_BAMLIST.out.bam2_text.map { [ it[0], it[1], it[2] ] }

        //
        // Run rMATS prep step
        //

        RMATS_PREP (
            gtf,
            ch_bamlist_cond1,
            ch_bamlist_cond2,
            ch_contrast_cond1,
            ch_contrast_cond2
        )

        //
        // Run rMATs post step
        //

        RMATS_POST (
            gtf,
            ch_bamlist_cond1,
            ch_bamlist_cond2,
            ch_contrast_cond1,
            ch_contrast_cond2,
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
