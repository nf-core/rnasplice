//
// tximport of the Salmon results of all samples
//

include { UNTAR            } from '../../../modules/nf-core/untar'
include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT as TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'
include { SUPPATPM         } from '../../../modules/local/suppatpm'

workflow TX2GENE_TXIMPORT {
    take:
    ch_salmon_results // channel: [ val(meta), path(results) ], a Salmon directory or a .tar.gz of one
    ch_gtf            // channel: path(gtf)
    gene_id           // string: GTF attribute of the gene ID
    gene_extra        // string: comma-separated GTF attributes to carry, e.g. the gene name

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

    // The Salmon directories of all samples: [ meta, [ results ] ]
    ch_merged_results = ch_extracted_results
        .map { _meta, results -> results }
        .collect()
        .map { results -> [[id: 'salmon.merged'], results] }

    //
    // MODULE: Transcript to gene table, with the GTF attribute that matches the Salmon transcript IDs
    //
    CUSTOM_TX2GENE(
        ch_gtf.map { gtf -> [[id: gtf.baseName], gtf] },
        ch_merged_results,
        'salmon',
        gene_id,
        gene_extra,
    )

    //
    // MODULE: Merge the quantification of all samples
    //
    TXIMPORT(ch_merged_results, CUSTOM_TX2GENE.out.tx2gene, 'salmon')

    //
    // MODULE: Transcript TPM table in the format SUPPA reads
    //
    SUPPATPM(TXIMPORT.out.tpm_transcript.map { _meta, tpm -> [[id: 'suppa_tpm'], tpm] })

    emit:
    salmon_results               = ch_extracted_results                    // channel: [ val(meta), path(results) ], the Salmon directory of each sample
    tx2gene                      = CUSTOM_TX2GENE.out.tx2gene              // channel: [ val(meta), path(tx2gene.tsv) ]
    tpm_transcript               = TXIMPORT.out.tpm_transcript             // channel: [ val(meta), path(transcript_tpm.tsv) ]
    counts_transcript            = TXIMPORT.out.counts_transcript          // channel: [ val(meta), path(transcript_counts.tsv) ]
    counts_transcript_scaled     = TXIMPORT.out.counts_transcript_scaled     // channel: [ val(meta), path(transcript_counts_scaled.tsv) ]
    counts_transcript_dtu_scaled = TXIMPORT.out.counts_transcript_dtu_scaled // channel: [ val(meta), path(transcript_counts_dtu_scaled.tsv) ]
    tpm_gene                     = TXIMPORT.out.tpm_gene                   // channel: [ val(meta), path(gene_tpm.tsv) ]
    counts_gene                  = TXIMPORT.out.counts_gene                // channel: [ val(meta), path(gene_counts.tsv) ]
    counts_gene_scaled           = TXIMPORT.out.counts_gene_scaled         // channel: [ val(meta), path(gene_counts_scaled.tsv) ]
    counts_gene_length_scaled    = TXIMPORT.out.counts_gene_length_scaled  // channel: [ val(meta), path(gene_counts_length_scaled.tsv) ]
    suppa_tpm                    = SUPPATPM.out.tpm.map { _meta, tpm -> tpm } // channel: path(suppa_tpm.txt)
}
