//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    format // value, get the format type of the samplesheet
    samplesheet_reformatted // channel for getting the re-formatted samplesheet as input

    main:

    switch(format) {
        case '[FASTQ]':
//            println("Input_check from 'fastq'");
            SAMPLESHEET_CHECK ( samplesheet_reformatted, format ).csv.splitCsv ( header:true, sep:',' ).map { create_fastq_channel(it) }.set { out }
            break;
        case '[TRANSCRIPTOME]':
//            println("Input_check from 'transcriptome_bam'");
            SAMPLESHEET_CHECK ( samplesheet_reformatted, format  ).csv.splitCsv ( header:true, sep:',' ).map { create_transcriptome_channel(it) }.set { out }
            break;
        case '[BAM]':
//            println("Input_check from 'bam'");
            SAMPLESHEET_CHECK ( samplesheet_reformatted, format ).csv.splitCsv ( header:true, sep:',' ).map { create_bam_channel(it) }.set { out }
            break;
        case '[SALMON]':
//            println("Input_check from 'salmon'");
            SAMPLESHEET_CHECK ( samplesheet_reformatted, format ).csv.splitCsv ( header:true, sep:',' ).map { create_salmon_channel(it) }.set { out }
            break;
        default: log.info("Doesn't work!!!!")
    }

    emit:
    out                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.strandedness = row.strandedness
    meta.condition = row.condition

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
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
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

// Function to get list of [ meta, bam ]
def create_transcriptome_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.condition = row.condition
    meta.bam = row.bam
    meta.transcriptome = row.transcriptome

    // add path(s) of the bam file(s) to the meta map
    def transcriptome_meta = []
    if (!file(row.transcriptome).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bam file does not exist!\n${row.transcriptome}"
    }
        transcriptome_meta = [ meta, [ file(row.transcriptome) ] ]
    return transcriptome_meta
}

// Function to get list of [ meta, salmon ]
def create_salmon_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.condition = row.condition
    meta.bam = row.salmon

    // add path(s) of the bam file(s) to the meta map
    def salmon_meta = []
    if (!file(row.salmon).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> salmon path does not exist!\n${row.salmon}"
    }
        salmon_meta = [ meta, [ file(row.salmon) ] ]
    return salmon_meta
}
