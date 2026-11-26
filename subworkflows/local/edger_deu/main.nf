//
// edgeR DEU subworkflow
//

include { SUBREAD_FLATTENGTF      } from '../../../modules/local/subread/flattengtf'
include { SUBREAD_FEATURECOUNTS   } from '../../../modules/nf-core/subread/featurecounts'
include { CSVTK_JOIN              } from '../../../modules/nf-core/csvtk/join'
include { EDGER_DIFFSPLICE        } from '../../../modules/nf-core/edger/diffsplice'
include { GENEANNOTATION_ANNOTATE } from '../../../modules/local/geneannotation/annotate'

workflow EDGER_DEU {
    take:
    gtf                // channel: path(gtf)
    ch_genome_bam      // channel: [ val(meta), path(bam) ]
    ch_contrastsheet   // channel: [ contrast:, treatment:, control:, blocking: (when set) ], the parsed contrasts
    ch_samplesheet     // channel: file(samplesheet), read for the values of the blocking variables
    ch_gene_annotation // channel: [ val(meta), path(tsv) ], the gene annotation table, as a value channel

    main:

    //
    // MODULE: SUBREAD_FLATTENGTF
    //

    SUBREAD_FLATTENGTF(gtf.map { file -> [[id: file.baseName], file] })

    //
    // MODULE: SUBREAD_FEATURECOUNTS
    //

    // Counted one sample at a time, so that each sample is counted with its own
    // strandedness and pairing
    ch_feature_counts = ch_genome_bam.combine(SUBREAD_FLATTENGTF.out.saf.map { _meta, saf -> saf })

    SUBREAD_FEATURECOUNTS(ch_feature_counts)

    //
    // MODULE: CSVTK_JOIN
    //

    CSVTK_JOIN(
        SUBREAD_FEATURECOUNTS.out.counts
            .map { _meta, counts -> counts }
            .toSortedList { a, b -> a.name <=> b.name }
            .map { counts -> [[id: 'edger_exon'], counts] }
    )

    // featureCounts names each count column after its BAM file, so the sample sheet the
    // model reads maps those names to the conditions and the blocking variables of all
    // contrasts. A sample with several sequencing runs takes its first row
    ch_blocking_columns = ch_contrastsheet
        .flatMap { row -> row.blocking ? row.blocking.tokenize(';') : [] }
        .unique()
        .toSortedList()

    ch_sample_rows = ch_samplesheet
        .splitCsv(header: true)
        .unique { row -> row.sample }
        .toList()

    ch_edger_samplesheet = ch_genome_bam
        .toList()
        .combine(ch_blocking_columns)
        .combine(ch_sample_rows)
        .map { bams, blocking_columns, rows ->
            def rows_by_sample = rows.collectEntries { row -> [ row.sample, row ] }
            def lines = [ (['sample', 'condition'] + blocking_columns).join(',') ]
            bams.sort { bam_a, bam_b -> bam_a[1].name <=> bam_b[1].name }.each { meta, bam ->
                lines << ([ bam.name, meta.condition ] + blocking_columns.collect { column -> rows_by_sample[meta.id][column] }).join(',')
            }
            lines.join('\n') + '\n'
        }
        .collectFile(name: 'edger_samplesheet.csv')
        .map { samplesheet -> [[id: 'edger_exon'], samplesheet] }

    //
    // MODULE: EDGER_DIFFSPLICE
    //

    ch_diffsplice_input = ch_contrastsheet
        .combine(ch_edger_samplesheet.join(CSVTK_JOIN.out.out_file))
        .multiMap { row, _meta, samplesheet, counts ->
            def meta = row.blocking ? [id: row.contrast, blocking: row.blocking] : [id: row.contrast]
            contrast: [meta, 'condition', row.control, row.treatment, '', '']
            counts: [meta, samplesheet, counts]
        }

    EDGER_DIFFSPLICE(
        ch_diffsplice_input.contrast,
        ch_diffsplice_input.counts,
    )

    //
    // MODULE: Copies of the edgeR tables with the gene names
    //

    GENEANNOTATION_ANNOTATE(
        EDGER_DIFFSPLICE.out.results.mix(EDGER_DIFFSPLICE.out.simes, EDGER_DIFFSPLICE.out.gene),
        ch_gene_annotation,
        'Geneid',
    )

    emit:
    featureCounts_counts  = SUBREAD_FEATURECOUNTS.out.counts   // channel: [ val(meta), path(counts) ]
    featureCounts_summary = SUBREAD_FEATURECOUNTS.out.summary  // channel: [ val(meta), path(summary) ]
    counts                = CSVTK_JOIN.out.out_file            // channel: [ val(meta), path(tsv) ], the counts of all samples
    results               = EDGER_DIFFSPLICE.out.results       // channel: [ val(meta), path(tsv) ], feature-level results per contrast
    simes                 = EDGER_DIFFSPLICE.out.simes         // channel: [ val(meta), path(tsv) ], gene-level Simes results per contrast
    gene                  = EDGER_DIFFSPLICE.out.gene          // channel: [ val(meta), path(tsv) ], gene-level test results per contrast
    plots                 = EDGER_DIFFSPLICE.out.plots         // channel: [ val(meta), path(pdf) ]
    bcv_plot              = EDGER_DIFFSPLICE.out.bcv_plot      // channel: [ val(meta), path(png) ]
    rdata                 = EDGER_DIFFSPLICE.out.rdata         // channel: [ val(meta), path(rds) ], the DGEGLM object
    lrt                   = EDGER_DIFFSPLICE.out.lrt           // channel: [ val(meta), path(rds) ], the DGELRT object of diffSpliceDGE
    annotated             = GENEANNOTATION_ANNOTATE.out.annotated // channel: [ val(meta), path(tsv) ], one table per item
}
