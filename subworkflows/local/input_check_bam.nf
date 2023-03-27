//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_bam_channel(it) }
        .set { bam }

    emit:
    bam                                       // channel: [ val(meta), [ bam ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, bam ] ]
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
//    meta.single_end = row.single_end.toBoolean()
//    meta.strandedness = row.strandedness
    meta.condition = row.condition
    meta.bam = row.bam

    // add path(s) of the bam file(s) to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bam file does not exist!\n${row.bam}"
    }
        bam_meta = [ meta, [ file(row.bam) ] ]
    return bam_meta
}

