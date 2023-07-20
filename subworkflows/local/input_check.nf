//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    source      // value: params.source

    main:

    switch(source) {
        case 'fastq':
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .set { reads }
            break;
        case 'genome_bam':
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_genome_bam_channel(it) }
            .set { reads }
            break;
        case 'transcriptome_bam':
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_transcriptome_bam_channel(it) }
            .set { reads }
            break;
        case 'salmon_results':
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_salmon_results_channel(it) }
            .set { reads }
            break;
    }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    meta.strandedness = row.strandedness
    meta.condition    = row.condition

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

// Function to get list of [ meta, bam ]
def create_genome_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.genome_bam = row.genome_bam
    meta.condition  = row.condition

    // add path(s) of the bam file(s) to the meta map
    def genome_bam_meta = []
    if (!file(row.genome_bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> genome bam file does not exist!\n${row.genome_bam}"
    }
    genome_bam_meta = [ meta, [ file(row.genome_bam) ] ]
    return genome_bam_meta
}

// Function to get list of [ meta, txbam ]
def create_transcriptome_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                = row.sample
    meta.transcriptome_bam = row.transcriptome_bam
    meta.condition         = row.condition

    // add path(s) of the bam and txbam file(s) to the meta map
    def transcriptome_bam_meta = []
    if (!file(row.transcriptome_bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> transcriptome bam file does not exist!\n${row.transcriptome_bam}"
    }
    transcriptome_bam_meta = [ meta, [ file(row.transcriptome_bam) ] ]
    return transcriptome_bam_meta
}

// Function to get list of [ meta, salmon ]
def create_salmon_results_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.sample
    meta.salmon_results = row.salmon_results
    meta.condition      = row.condition

    // add path(s) of the salmon file(s) to the meta map
    def salmon_results_meta = []
    if (!file(row.salmon_results).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> salmon results directory does not exist!\n${row.salmon_results}"
    }
    salmon_results_meta = [ meta, [ file(row.salmon_results) ] ]
    return salmon_results_meta
}
