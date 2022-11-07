//
// Read QC, and trimming
//

include { FASTQC           } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE       } from '../../modules/nf-core/trimgalore/main'

workflow FASTQC_TRIMGALORE {

    take:

    reads            // channel: [ val(meta), [ reads ] ]
    skip_fastqc      // boolean: true/false
    skip_trimming    // boolean: true/false

    main:

    // define empty versions channel
    ch_versions = Channel.empty()

    // define fastqc output channels
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()

    if (!skip_fastqc) {

        // Run FASTQC module on reads channel [ val(meta), [ reads ] ]
        FASTQC ( reads )

        // take FASTQC output and set as separate channels
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip

        // Add software versions to versions channel
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    // define trimgalore output channels
    trim_reads    = Channel.empty()
    trim_unpaired = Channel.empty()
    trim_html     = Channel.empty()
    trim_zip      = Channel.empty()
    trim_log      = Channel.empty()

    if (!skip_trimming) {

        // Run TRIMGALORE module on reads channel [ val(meta), [ reads ] ]
        TRIMGALORE ( reads )

        // take TRIMGALORE output and set as separate channels
        trim_reads    = TRIMGALORE.out.reads
        trim_unpaired = TRIMGALORE.out.unpaired
        trim_html     = TRIMGALORE.out.html
        trim_zip      = TRIMGALORE.out.zip
        trim_log      = TRIMGALORE.out.log

        // Add software versions to versions channel
        ch_versions   = ch_versions.mix(TRIMGALORE.out.versions.first())
    }

    // Define output
    emit:

    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    trim_unpaired      // channel: [ val(meta), [ reads ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
