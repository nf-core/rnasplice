//
// Visualise splicing events with MISO sashimi plots
//

include { GFFREAD as GFFREAD_GTF2GFF3 } from '../../../modules/nf-core/gffread/main'
include { MISOPY_INDEX                } from '../../../modules/nf-core/misopy/index/main'
include { MISOPY_RUN                  } from '../../../modules/local/misopy/run/main'
include { MISOPY_SASHIMIPLOT          } from '../../../modules/local/misopy/sashimiplot/main'
include { MISOPY_SETTINGS             } from '../../../modules/local/misopysettings/main'

workflow VISUALISE_MISO {

    take:
    gtf             // channel: [ path(gtf) ]
    ch_genome_bam   // channel: [ val(meta), path(bam) ]
    ch_genome_bai   // channel: [ val(meta), path(bai) ]
    fig_width       // val: figure width in inches
    fig_height      // val: figure height in inches
    miso_genes      // val: comma separated list of events to plot
    miso_genes_file // val: path to a file listing events to plot, or false

    main:

    //
    // MODULE: GFFREAD_GTF2GFF3
    //

    GFFREAD_GTF2GFF3(gtf.map { annotation -> [[id: annotation.baseName], annotation] }, [])

    //
    // MODULE: MISOPY_INDEX
    //

    MISOPY_INDEX(GFFREAD_GTF2GFF3.out.gffread_gff)

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
    // MODULE: MISOPY_SETTINGS
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
    // MODULE: MISOPY_SASHIMIPLOT
    //

    def miso_genes_list = miso_genes ? miso_genes.split(',').collect { it -> it.trim() } : []
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

    // sashimi_plot fails on a gene the index does not have, so requested genes are matched to
    // the index file names, which MISO takes from the GFF3 gene IDs, ignoring version suffixes
    ch_miso_events = ch_miso_index
        .combine(ch_miso_genes.collect().map { genes -> [genes] })
        .flatMap { _meta, index, genes ->
            def index_genes = index.listFiles()
                .findAll { chrom -> chrom.isDirectory() }
                .collectMany { chrom -> chrom.listFiles().findAll { gene -> gene.extension == 'pickle' }.collect { gene -> gene.baseName } }
            genes.collect { gene ->
                def unversioned = gene.replaceFirst(/\.\d+$/, '')
                def event = index_genes.find { id -> id == gene } ?: index_genes.find { id -> id.replaceFirst(/\.\d+$/, '') == unversioned }
                if (!event) {
                    log.warn("No sashimi plot for '${gene}': the GTF has no gene with that ID, with or without a version suffix.")
                }
                else if (event != gene && gene != unversioned) {
                    log.warn("Plotting '${event}' for '${gene}': the GTF has that gene under a different version.")
                }
                return event
            }.findAll { event -> event }.unique()
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
        .combine(ch_miso_events)

    MISOPY_SASHIMIPLOT(
        ch_sashimiplot_input,
        ch_miso_index,
    )

    emit:
    gff3          = GFFREAD_GTF2GFF3.out.gffread_gff   // channel: [ val(meta), path(*.gff3) ]
    miso_index    = ch_miso_index                      // channel: [ val(meta), path(index) ]
    miso_data     = ch_miso_data                       // channel: [ val(meta), path(miso_data/*) ]
    miso_settings = MISOPY_SETTINGS.out.miso_settings  // channel: [ val(meta), path(miso_settings.txt) ]
    miso_sashimi  = MISOPY_SASHIMIPLOT.out.sashimi_plot // channel: [ val(meta), path(*.pdf) ]
}
