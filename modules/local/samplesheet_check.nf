process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path samplesheet
    val source

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnasplice/bin/

    switch (source) {
        case 'fastq':
            """
            check_samplesheet_fastq.py $samplesheet samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
        case 'genome_bam':
            """
            check_samplesheet_genome_bam.py $samplesheet samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
        case 'transcriptome_bam':
            """
            check_samplesheet_transcriptome_bam.py $samplesheet samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
        case 'salmon_results':
            """
            check_samplesheet_salmon_results.py $samplesheet samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
    }

}
