//
// Infer the strandedness of the samples set to 'auto' in the samplesheet
//
// This follows nf-core/rnaseq 3.22.2, where the reads of an 'auto' sample are sub-sampled
// and pseudo-aligned with Salmon in its own library type detection mode, and the library
// type counts are called forward, reverse or unstranded:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
//

include { FASTQ_SUBSAMPLE_FQ_SALMON     } from '../../nf-core/fastq_subsample_fq_salmon'
include { getSalmonInferredStrandedness } from '../utils_nfcore_rnasplice_pipeline'

workflow INFER_STRANDEDNESS {
    take:
    ch_reads             // channel: [ val(meta), [ reads ] ], meta.strandedness may be 'auto'
    ch_fasta             // channel: path(genome.fasta)
    ch_transcript_fasta  // channel: path(transcript.fasta)
    ch_gtf               // channel: path(genome.gtf)
    ch_salmon_index      // channel: path(salmon/index/)
    make_salmon_index    // boolean: build a Salmon index, as the pipeline has not been given or built one
    stranded_threshold   // float: fraction of stranded reads a strand must reach to be called forward or reverse
    unstranded_threshold // float: difference between the forward and reverse fractions below which a sample is called unstranded

    main:

    ch_reads
        .branch { meta, reads ->
            auto_strand: meta.strandedness == 'auto'
            return [meta, reads]
            known_strand: meta.strandedness != 'auto'
            return [meta, reads]
        }
        .set { ch_strand_reads }

    // The genome is only handed on once there is an 'auto' sample, so that a Salmon index
    // is never built for a run that does not need one
    ch_genome_fasta = ch_fasta
        .combine(ch_strand_reads.auto_strand)
        .map { items -> items.first() }
        .first()

    FASTQ_SUBSAMPLE_FQ_SALMON(
        ch_strand_reads.auto_strand,
        ch_genome_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        make_salmon_index,
    )

    // An undetermined sample is analysed as unstranded, as nf-core/rnaseq does. Its
    // fractions stay in the meta map so that the RSeQC strand check can report them
    FASTQ_SUBSAMPLE_FQ_SALMON.out.lib_format_counts
        .join(ch_strand_reads.auto_strand)
        .map { meta, json, reads ->
            def salmon_strand_analysis = getSalmonInferredStrandedness(json, stranded_threshold, unstranded_threshold)
            def strandedness = salmon_strand_analysis.inferred_strandedness
            if (strandedness == 'undetermined') {
                log.warn("Strandedness of sample '${meta.id}' could not be inferred with Salmon, so it is analysed as 'unstranded'. Set it in the samplesheet if that is wrong.")
                strandedness = 'unstranded'
            }
            return [meta + [strandedness: strandedness, salmon_strand_analysis: salmon_strand_analysis], reads]
        }
        .mix(ch_strand_reads.known_strand)
        .set { ch_strand_inferred_reads }

    emit:
    reads             = ch_strand_inferred_reads                            // channel: [ val(meta), [ reads ] ]
    lib_format_counts = FASTQ_SUBSAMPLE_FQ_SALMON.out.lib_format_counts     // channel: [ val(meta), path(lib_format_counts.json) ]
}
