//
// Visualise miso subworkflow
//

include { GTF_2_GFF3    } from '../../modules/local/gtf_2_gff3'
include { MISO_INDEX    } from '../../modules/local/miso_index'
include { MISO_RUN      } from '../../modules/local/miso_run'
include { MISO_SETTINGS } from '../../modules/local/miso_settings'
include { MISO_SASHIMI  } from '../../modules/local/miso_sashimi'


workflow VISUALISE_MISO {

    take:

    gtf                // path gtf
    ch_genome_bam      // channel: [ val(meta), path(bams) ]
    ch_genome_bai      // channel: [ val(meta), path(bais) ]
    miso_read_len      // 75
    fig_width          // 7
    fig_height         // 5
    miso_genes         // params.miso_genes
    miso_genes_file    // params.miso_genes_file

    main:

    ch_versions = Channel.empty()

   //
   // MODULE: gtf_2_gff3
   //

    GTF_2_GFF3 (
        gtf
    )

    ch_versions = ch_versions.mix(GTF_2_GFF3.out.versions)

   //
   // MODULE: DEXSeq Annotation
   //

    def index_prefix = "index"

    MISO_INDEX (
        GTF_2_GFF3.out.gff3,
        index_prefix
    )

    ch_versions = ch_versions.mix(MISO_INDEX.out.versions)

    //
    // MODULE: MISO_RUN
    //

    ch_bam_join = ch_genome_bam.join(ch_genome_bai)
    ch_miso_index =  MISO_INDEX.out.miso_index.collect()

    MISO_RUN (
        ch_miso_index,
        ch_bam_join,
        miso_read_len
    )

    ch_versions = ch_versions.mix(MISO_RUN.out.versions)

    //
    // MODULE: MISO_SETTINGS
    //

    ch_bams = ch_genome_bam.collect({it[1]})
    ch_miso_run = MISO_RUN.out.miso.map{it[1]}.collect()

    MISO_SETTINGS (
        ch_miso_run,
        ch_bams,
        fig_width,
        fig_height
    )

    ch_versions = ch_versions.mix(MISO_SETTINGS.out.versions)

    //
    // MODULE: MISO_SASHIMI
    //

    def miso_genes_list = miso_genes ? miso_genes.split(',').collect{ it.trim() } : [""]
    ch_miso_genes_list = Channel.fromList( miso_genes_list )

    if (miso_genes_file && miso_genes) {
        ch_miso_genes_file = Channel.fromPath(miso_genes_file)
            .splitCsv()
        ch_miso_genes_list.concat( ch_miso_genes_file )
        .set{ ch_miso_genes }
    } else if (miso_genes_file) {
        ch_miso_genes_file = Channel.fromPath(miso_genes_file)
            .splitCsv()
            .set{ ch_miso_genes }
        } else {
        ch_miso_genes = ch_miso_genes_list
    }
    ch_miso_input = MISO_SETTINGS.out.miso_settings.combine(ch_miso_genes)

    ch_bam_bai = ch_bam_join
        .map { [it[1], it[2] ] }
        .collect()

    MISO_SASHIMI (
        ch_miso_index,
        ch_miso_input,
        ch_bam_bai,
        ch_miso_run
    )

    ch_versions = ch_versions.mix(MISO_SASHIMI.out.versions)

    emit:

    gff3                        = GTF_2_GFF3.out.gff3               // path *.gff3
    miso_index                  = ch_miso_index                     // channel: [ ch_miso_index ]
    miso_data                   = ch_miso_run                       // channel: [ ch_miso_run ]
    miso_settings               = MISO_SETTINGS.out.miso_settings   // path miso_setting.txt
    miso_sashimi                = MISO_SASHIMI.out.sashimi          // path sashimi/*.pdf

    versions                    = ch_versions                  // channel: [ versions.yml ]
}
