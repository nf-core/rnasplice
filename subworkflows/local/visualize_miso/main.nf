//
// Visualise miso subworkflow
//

include { GTF_2_GFF3         } from '../../../modules/local/gtf_2_gff3'
include { MISOPY_INDEX       } from '../../../modules/nf-core/misopy/index'
include { MISOPY_RUN         } from '../../../modules/local/misopy/run'
include { MISOPY_SASHIMIPLOT } from '../../../modules/local/misopy/sashimiplot'
include { MISOPY_SETTINGS    } from '../../../modules/local/misopysettings'


workflow VISUALISE_MISO {
    take:
    gtf // path gtf
    ch_genome_bam // channel: [ val(meta), path(bam) ]
    ch_genome_bai // channel: [ val(meta), path(bai) ]
    fig_width
    fig_height
    miso_genes // params.miso_genes
    miso_genes_file // params.miso_genes_file

    main:

    //
    // MODULE: gtf_2_gff3
    //

    GTF_2_GFF3(gtf)

    //
    // MODULE: DEXSeq Annotation
    //

    MISOPY_INDEX(GTF_2_GFF3.out.gff3.map { gff -> [[:], gff] })

    //
    // MODULE: MISOPY_RUN
    //

    ch_bam_bai = ch_genome_bam.join(ch_genome_bai)
    ch_miso_index = MISOPY_INDEX.out.miso_index

    MISOPY_RUN(
        ch_bam_bai,
        ch_miso_index,
    )

    ch_miso_data = MISOPY_RUN.out.miso
        .collect { _meta, miso_data -> miso_data }
        .map { miso_data -> [[id: 'miso'], miso_data] }

    //
    // MODULE: MISO_SETTINGS
    //

    ch_miso_settings_input = ch_genome_bam
        .collect { _meta, bam -> bam }
        .map { bams -> [[id: 'miso'], bams] }
        .join(ch_miso_data, by: 0)

    ch_miso_settings = MISOPY_SETTINGS(
        ch_miso_settings_input,
        fig_width,
        fig_height,
    ).miso_settings

    //
    // MODULE: MISO_SASHIMI
    //

    def miso_genes_list = miso_genes ? miso_genes.split(',').collect { it -> it.trim() } : [""]
    ch_miso_genes_list = channel.fromList(miso_genes_list)

    if (miso_genes_file && miso_genes) {
        ch_miso_genes_file = channel.fromPath(miso_genes_file)
            .splitCsv()
        ch_miso_genes = ch_miso_genes_list.concat(ch_miso_genes_file)
    }
    else if (miso_genes_file) {
        ch_miso_genes = channel.fromPath(miso_genes_file)
            .splitCsv()
    }
    else {
        ch_miso_genes = ch_miso_genes_list
    }

    ch_bams = ch_genome_bam
        .collect { _meta, bam -> bam }
        .map { bams -> [[id: 'miso'], bams] }

    ch_bais = ch_genome_bai
        .collect { _meta, bai -> bai }
        .map { bais -> [[id: 'miso'], bais] }

    ch_sashimiplot_input = ch_bams
        .join(ch_bais)
        .join(ch_miso_data, by: 0)
        .join(ch_miso_settings, by: 0)
        .combine(ch_miso_genes)

    MISOPY_SASHIMIPLOT(
        ch_sashimiplot_input,
        ch_miso_index,
    )

    emit:
    gff3          = GTF_2_GFF3.out.gff3 // path *.gff3
    miso_index    = ch_miso_index // channel: [ ch_miso_index ]
    miso_data     = ch_miso_data // channel: [ ch_miso_data ]
    miso_settings = MISOPY_SETTINGS.out.miso_settings // path miso_setting.txt
    miso_sashimi  = MISOPY_SASHIMIPLOT.out.sashimi_plot // path sashimi/*.pdf
}
