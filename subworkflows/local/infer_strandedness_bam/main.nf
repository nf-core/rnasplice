//
// Infer the strandedness of the genome BAM samples set to 'auto' in the samplesheet
//
// There are no reads to pseudo-align with Salmon here, so the strandedness is called from
// the RSeQC infer_experiment.py report of the BAM, with the same thresholds and the same
// parser nf-core/rnaseq uses for its RSeQC strand check:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/local/utils_nfcore_rnaseq_pipeline/main.nf
//

include { RSEQC_INFEREXPERIMENT          } from '../../../modules/nf-core/rseqc/inferexperiment'
include { getInferexperimentStrandedness } from '../infer_strandedness'

workflow INFER_STRANDEDNESS_BAM {
    take:
    ch_bam_bai           // channel: [ val(meta), path(bam), path(bai) ], meta.strandedness may be 'auto'
    ch_bed               // channel: path(genes.bed), the BED12 annotation RSeQC needs
    stranded_threshold   // float: fraction of stranded reads a strand must reach to be called forward or reverse
    unstranded_threshold // float: difference between the forward and reverse fractions below which a sample is called unstranded

    main:

    ch_bam_bai
        .branch { meta, bam, bai ->
            auto_strand: meta.strandedness == 'auto'
            return [meta, bam, bai]
            known_strand: meta.strandedness != 'auto'
            return [meta, bam, bai]
        }
        .set { ch_strand_bam }

    RSEQC_INFEREXPERIMENT(
        ch_strand_bam.auto_strand,
        ch_bed,
    )

    // An undetermined sample is analysed as unstranded, as for the reads. Its fractions stay
    // in the meta map so that the strand check can report them without running RSeQC again
    RSEQC_INFEREXPERIMENT.out.txt
        .join(ch_strand_bam.auto_strand)
        .map { meta, txt, bam, bai ->
            def rseqc_strand_analysis = getInferexperimentStrandedness(txt, stranded_threshold, unstranded_threshold)
            def strandedness = rseqc_strand_analysis.inferred_strandedness
            if (strandedness == 'undetermined') {
                log.warn("Strandedness of sample '${meta.id}' could not be inferred with RSeQC, so it is analysed as 'unstranded'. Set it in the samplesheet if that is wrong.")
                strandedness = 'unstranded'
            }
            return [meta + [strandedness: strandedness, rseqc_strand_analysis: rseqc_strand_analysis], bam, bai, txt]
        }
        .set { ch_strand_inferred }

    emit:
    bam_bai = ch_strand_inferred.map { meta, bam, bai, _txt -> [meta, bam, bai] }.mix(ch_strand_bam.known_strand) // channel: [ val(meta), path(bam), path(bai) ]
    txt     = ch_strand_inferred.map { meta, _bam, _bai, txt -> [meta, txt] }                                       // channel: [ val(meta), path(infer_experiment.txt) ]
}
