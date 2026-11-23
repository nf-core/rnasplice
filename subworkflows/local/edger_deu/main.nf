//
// edgeR DEU subworkflow
//

include { SUBREAD_FLATTENGTF } from '../../../modules/local/subread/flattengtf'
include { SUBREAD_FEATURECOUNTS } from '../../../modules/nf-core/subread/featurecounts'
include { EDGER_EXON } from '../../../modules/local/edger/exon'
include { GENEANNOTATION_ANNOTATE } from '../../../modules/local/geneannotation/annotate'

workflow EDGER_DEU {
    take:
    gtf // path: gtf
    ch_genome_bam // channel: [ val(meta), path(bams) ]
    ch_samplesheet // channel.fromPath(params.input)
    ch_contrastsheet // channel.fromPath(params.contrasts)
    ch_gene_annotation // channel: [ val(meta), path(tsv) ], the gene annotation table, as a value channel

    main:

    //
    // MODULE: SUBREAD_FLATTENGTF
    //

    SUBREAD_FLATTENGTF(gtf.map { file -> [[id: file.baseName], file] })

    //
    // MODULE: SUBREAD_FEATURECOUNTS
    //

    ch_feature_counts = ch_genome_bam.combine(SUBREAD_FLATTENGTF.out.saf.map { _meta, saf -> saf })

    SUBREAD_FEATURECOUNTS(ch_feature_counts)

    //
    // MODULE: EDGER_EXON
    //

    ch_feature_counts_collected = SUBREAD_FEATURECOUNTS.out.counts
        .map { _meta, counts -> counts }
        .collect()
        .map { counts -> [[id: 'edger_exon'], counts] }

    EDGER_EXON(
        ch_feature_counts_collected,
        ch_samplesheet.map { samplesheet -> [[id: samplesheet.baseName], samplesheet] },
        ch_contrastsheet.map { contrastsheet -> [[id: contrastsheet.baseName], contrastsheet] },
    )

    //
    // MODULE: Copies of the edgeR tables with the gene names
    //

    // EDGER_EXON writes the tables of every contrast in one task, which are annotated
    // one by one
    GENEANNOTATION_ANNOTATE(
        EDGER_EXON.out.edger_exon_csv.transpose(),
        ch_gene_annotation,
        'Geneid',
    )

    emit:
    featureCounts_counts = SUBREAD_FEATURECOUNTS.out.counts // channel: [ val(meta), path(counts) ]
    featureCounts_summary = SUBREAD_FEATURECOUNTS.out.summary // channel: [ val(meta), path(summary) ]
    edger_exon_dge = EDGER_EXON.out.edger_exon_dge // channel: [ val(meta), path(rds) ]
    edger_exon_glm = EDGER_EXON.out.edger_exon_glm // channel: [ val(meta), path(rds) ]
    edger_exon_lrt = EDGER_EXON.out.edger_exon_lrt // channel: [ val(meta), path(rds) ]
    edger_exon_csv = EDGER_EXON.out.edger_exon_csv // channel: [ val(meta), path(csv) ]
    edger_exon_pdf = EDGER_EXON.out.edger_exon_pdf // channel: [ val(meta), path(pdf) ]
    edger_exon_annotated_csv = GENEANNOTATION_ANNOTATE.out.annotated // channel: [ val(meta), path(csv) ], one table per item
}
