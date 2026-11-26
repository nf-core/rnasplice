//
// tximport of the Salmon results of all samples
//

include { UNTAR                               } from '../../../modules/nf-core/untar'
include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../../nf-core/quant_tximport_summarizedexperiment'
include { SUPPATPM                            } from '../../../modules/local/suppatpm'

workflow TX2GENE_TXIMPORT {
    take:
    ch_samplesheet    // channel: [ val(meta), path(samplesheet) ], for the colData of the SummarizedExperiment objects
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

    //
    // SUBWORKFLOW: Transcript to gene table, tximport and SummarizedExperiment objects of all samples
    //
    QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        ch_samplesheet,
        ch_extracted_results,
        ch_gtf,
        gene_id,
        gene_extra,
        'salmon',
        false,
    )
    ch_txi = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT.out

    //
    // MODULE: Transcript TPM table in the format SUPPA reads
    //
    SUPPATPM(ch_txi.tpm_transcript.map { _meta, tpm -> [[id: 'suppa_tpm'], tpm] })

    emit:
    salmon_results               = ch_extracted_results                // channel: [ val(meta), path(results) ], the Salmon directory of each sample
    counts_transcript_scaled     = ch_txi.counts_transcript_scaled     // channel: [ val(meta), path(transcript_counts_scaled.tsv) ]
    counts_transcript_dtu_scaled = ch_txi.counts_transcript_dtu_scaled // channel: [ val(meta), path(transcript_counts_dtu_scaled.tsv) ]
    suppa_tpm                    = SUPPATPM.out.tpm.map { _meta, tpm -> tpm } // channel: path(suppa_tpm.txt)
}
