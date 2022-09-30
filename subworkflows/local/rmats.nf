//
// Uncompress and prepare reference genome files
//

include { CREATE_BAMLISTS                   } from '../../modules/local/create_bamlists.nf'
include { RMATS_PREP                        } from '../../modules/local/rmats_prep.nf'
include { RMATS_POST                        } from '../../modules/local/rmats_post.nf'

workflow RMATS {

    take:

    ch_genome_bam_conditions             // channel: [ [condition1, [condition1_metas], [condition1_bams]], [condition2, [condition2_metas], [condition2_bams]]]
    gtf                                  // channel: /path/to/genome.gtf

    main:

    ch_versions = Channel.empty()

    //
    // Split channel by meta condition
    //

    ch_genome_bam_conditions
      .multiMap { it -> 
        bam1: it[0]
        bam2: it[1]
      }
      .set { ch_genome_bam_conditions }

    //
    // Create input bam file
    //

    CREATE_BAMLISTS ( 
      ch_genome_bam_conditions.bam1,
      ch_genome_bam_conditions.bam2
    )
    

    //
    // Run rMATS prep step
    //

    RMATS_PREP ( 
      gtf,
      CREATE_BAMLISTS.out.bam_group1,
      CREATE_BAMLISTS.out.bam_group2,
      ch_genome_bam_conditions.bam1,
      ch_genome_bam_conditions.bam2 
    )

    //
    // Run rMATs post step
    //
    
    RMATS_POST ( 
      gtf,
      CREATE_BAMLISTS.out.bam_group1,
      CREATE_BAMLISTS.out.bam_group2,
      ch_genome_bam_conditions.bam1,
      ch_genome_bam_conditions.bam2,
      RMATS_PREP.out.rmats_temp
    )
    
    ch_versions    = ch_versions.mix(RMATS_POST.out.versions)
    
    // Define output

    emit:

    rmats_prep       = RMATS_PREP.out.rmats_temp         //    path: rmats_temp/*
    rmats_prep_log   = RMATS_PREP.out.log                //    path: rmats_prep.log
    rmats_post       = RMATS_POST.out.rmats_post         //    path: rmats_post/*
    rmats_post_log   = RMATS_POST.out.log                //    path: rmats_post.log

    versions = ch_versions.ifEmpty(null)                 //    channel: [ versions.yml ]
}

