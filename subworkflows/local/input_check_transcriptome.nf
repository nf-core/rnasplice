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
        .map { create_transcriptome_channel(it) }
        .set { transcriptome }

    emit:
    transcriptome                                       // channel: [ val(meta), [ transcriptome ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}


// Function to get list of [ meta, bam ] ]
def create_transcriptome_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
//    meta.single_end = row.single_end.toBoolean()
//    meta.strandedness = row.strandedness
    meta.condition = row.condition
    meta.bam = row.bam
    meta.transcriptome = row.transcriptome

    // add path(s) of the bam file(s) to the meta map
    def transcriptome_meta = []
    if (file(row.transcriptome).exists()) {
        transcriptome_meta = [ meta, [ file(row.transcriptome) ] ]
    return transcriptome_meta }
}
