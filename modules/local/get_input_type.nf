workflow GET_INPUT_TYPE {
    take:
    samplesheet

    main:
    myFile = samplesheet
    allLines = myFile.readLines()
    header = allLines.get(0)
    fields = header.split(',')

    switch (fields) {
        case {it.contains("fastq_1")}:
            string = "FASTQ"
            break;
        case {it.contains("transcriptome")}:
            string = "TRANSCRIPTOME"
            break;
        case {it.contains("bam")}:
            string = "BAM"
            break;
        case {it.contains("salmon")}:
            string = "SALMON"
            break;
        default: log.info("The samplesheet doesn't have a correct format!!!!")
        }

    emit:
    string

    }
