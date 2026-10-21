//
// Infer the strandedness of the samples set to 'auto' in the samplesheet
//
// This follows nf-core/rnaseq 3.22.2, where the reads of an 'auto' sample are sub-sampled
// and pseudo-aligned with Salmon in its own library type detection mode, and the library
// type counts are called forward, reverse or unstranded:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
//

include { FASTQ_SUBSAMPLE_FQ_SALMON } from '../../nf-core/fastq_subsample_fq_salmon'

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

//
// Call a library forward, reverse or unstranded from the reads assigned to each strand.
// Copied from nf-core/rnaseq 3.22.2:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
//
def calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def totalFragments = forwardFragments + reverseFragments + unstrandedFragments
    def totalStrandedFragments = forwardFragments + reverseFragments

    def strandedness = 'undetermined'
    if (totalStrandedFragments > 0) {
        def forwardProportion = forwardFragments / (totalStrandedFragments as double)
        def reverseProportion = reverseFragments / (totalStrandedFragments as double)
        def proportionDifference = Math.abs(forwardProportion - reverseProportion)

        if (forwardProportion >= stranded_threshold) {
            strandedness = 'forward'
        }
        else if (reverseProportion >= stranded_threshold) {
            strandedness = 'reverse'
        }
        else if (proportionDifference <= unstranded_threshold) {
            strandedness = 'unstranded'
        }
    }

    // Without any fragment, as from a stub run, the percentages are 0 rather than NaN,
    // which cannot be written to JSON
    def scale = totalFragments > 0 ? 100 / (totalFragments as double) : 0

    return [
        inferred_strandedness: strandedness,
        forwardFragments: forwardFragments * scale,
        reverseFragments: reverseFragments * scale,
        unstrandedFragments: unstrandedFragments * scale,
    ]
}

//
// Parse the Salmon 'lib_format_counts.json' file to get the inferred strandedness.
// Copied from nf-core/rnaseq 3.22.2:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
//
def getSalmonInferredStrandedness(json_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    // The Salmon stub writes an empty file, which is read as no counts at all
    def json_text = json_file.text
    def libCounts = json_text.trim() ? new groovy.json.JsonSlurper().parseText(json_text) : [:]

    def forwardKeys = ['SF', 'ISF', 'MSF', 'OSF']
    def reverseKeys = ['SR', 'ISR', 'MSR', 'OSR']

    // Every fragment Salmon assigns has a strand, even in an unstranded library, so these
    // counts are in practice always zero: https://groups.google.com/g/sailfish-users/c/yxzBDv6NB6I
    def unstrandedKeys = ['IU', 'U', 'MU']

    def forwardFragments = forwardKeys.collect { key -> libCounts[key] ?: 0 }.sum()
    def reverseFragments = reverseKeys.collect { key -> libCounts[key] ?: 0 }.sum()
    def unstrandedFragments = unstrandedKeys.collect { key -> libCounts[key] ?: 0 }.sum()

    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

//
// Parse the RSeQC 'infer_experiment.py' report to get the inferred strandedness, for
// both the single end and the paired end patterns. Copied from nf-core/rnaseq 3.22.2:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/local/utils_nfcore_rnaseq_pipeline/main.nf
//
def getInferexperimentStrandedness(inferexperiment_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def forwardFragments = 0
    def reverseFragments = 0
    def unstrandedFragments = 0

    inferexperiment_file.eachLine { line ->
        def unstranded_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
        def se_sense_matcher = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
        def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
        def pe_sense_matcher = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
        def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/

        if (unstranded_matcher) unstrandedFragments = unstranded_matcher[0][1].toFloat() * 100
        if (se_sense_matcher) forwardFragments = se_sense_matcher[0][1].toFloat() * 100
        if (se_antisense_matcher) reverseFragments = se_antisense_matcher[0][1].toFloat() * 100
        if (pe_sense_matcher) forwardFragments = pe_sense_matcher[0][1].toFloat() * 100
        if (pe_antisense_matcher) reverseFragments = pe_antisense_matcher[0][1].toFloat() * 100
    }

    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}
