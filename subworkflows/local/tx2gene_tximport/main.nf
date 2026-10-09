//
// Transcript to gene map and tximport of the Salmon results of all samples
//

include { GFFREAD_TX2GENE } from '../../../modules/local/gffread/tx2gene'
include { TXIMETA_TXIMPORT as TXIMPORT } from '../../../modules/local/tximeta/tximport'
include { UNTAR } from '../../../modules/nf-core/untar'

workflow TX2GENE_TXIMPORT {
    take:
    ch_salmon_results // channel: [ val(meta), path(results) ], a Salmon directory or a .tar.gz of one
    ch_gtf // channel: path(gtf)

    main:

    ch_branched_results = ch_salmon_results.branch { _meta, results ->
        tar: results.name.endsWith('.tar.gz')
        dir: true
    }

    //
    // MODULE: Extract the tarballs
    //
    UNTAR(ch_branched_results.tar)

    ch_extracted_results = ch_branched_results.dir.mix(UNTAR.out.untar)

    //
    // MODULE: Transcript to gene map of the annotation
    //
    GFFREAD_TX2GENE(ch_gtf)

    // The Salmon directories of all samples: [ meta, [ results ] ]
    ch_merged_results = ch_extracted_results
        .map { _meta, results -> results }
        .collect()
        .map { results -> [[id: 'salmon.merged'], results] }

    //
    // MODULE: Merge the quantification of all samples
    //
    TXIMPORT(ch_merged_results, GFFREAD_TX2GENE.out.tx2gene)

    emit:
    salmon_results = ch_extracted_results // channel: [ val(meta), path(results) ], the Salmon directory of each sample
    tx2gene = GFFREAD_TX2GENE.out.tx2gene // channel: path(tx2gene.tsv)
    txi = TXIMPORT.out.txi // channel: path(txi.rds)
    txi_s = TXIMPORT.out.txi_s // channel: path(txi.s.rds)
    txi_ls = TXIMPORT.out.txi_ls // channel: path(txi.ls.rds)
    txi_dtu = TXIMPORT.out.txi_dtu // channel: path(txi.dtu.rds)
    gi = TXIMPORT.out.gi // channel: path(gi.rds)
    gi_s = TXIMPORT.out.gi_s // channel: path(gi.s.rds)
    gi_ls = TXIMPORT.out.gi_ls // channel: path(gi.ls.rds)
    tpm_gene = TXIMPORT.out.tpm_gene // channel: path(gene_tpm.tsv)
    counts_gene = TXIMPORT.out.counts_gene // channel: path(gene_counts.tsv)
    tpm_gene_scaled = TXIMPORT.out.tpm_gene_scaled // channel: path(gene_tpm_scaled.tsv)
    counts_gene_scaled = TXIMPORT.out.counts_gene_scaled // channel: path(gene_counts_scaled.tsv)
    tpm_gene_length_scaled = TXIMPORT.out.tpm_gene_length_scaled // channel: path(gene_tpm_length_scaled.tsv)
    counts_gene_length_scaled = TXIMPORT.out.counts_gene_length_scaled // channel: path(gene_counts_length_scaled.tsv)
    tpm_transcript = TXIMPORT.out.tpm_transcript // channel: path(transcript_tpm.tsv)
    counts_transcript = TXIMPORT.out.counts_transcript // channel: path(transcript_counts.tsv)
    tpm_transcript_scaled = TXIMPORT.out.tpm_transcript_scaled // channel: path(transcript_tpm_scaled.tsv)
    counts_transcript_scaled = TXIMPORT.out.counts_transcript_scaled // channel: path(transcript_counts_scaled.tsv)
    tpm_transcript_length_scaled = TXIMPORT.out.tpm_transcript_length_scaled // channel: path(transcript_tpm_length_scaled.tsv)
    counts_transcript_length_scaled = TXIMPORT.out.counts_transcript_length_scaled // channel: path(transcript_counts_length_scaled.tsv)
    tpm_transcript_dtu_scaled = TXIMPORT.out.tpm_transcript_dtu_scaled // channel: path(transcript_tpm_dtu_scaled.tsv)
    counts_transcript_dtu_scaled = TXIMPORT.out.counts_transcript_dtu_scaled // channel: path(transcript_counts_dtu_scaled.tsv)
    tximport_tx2gene = TXIMPORT.out.tximport_tx2gene // channel: path(tximport.tx2gene.tsv)
    suppa_tpm = TXIMPORT.out.suppa_tpm // channel: path(suppa_tpm.txt)
}
