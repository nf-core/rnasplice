process SAMPLESHEET_REFORMAT {
    tag "$params.input"
    label 'process_single'

    conda "anaconda::pandas=1.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/wslh-bioinformatics/pandas:1.5.0-wslh-signed' }"

    input:
    path input

    output:
    path ("*.fastq.csv"), emit: fastq, optional: true
    path ("*.bam.csv"), emit: bam, optional: true
    path ("*.transcriptome.csv"), emit: transcriptome, optional: true
    path ("*.salmon.csv"), emit: salmon, optional: true
    path "versions.yml",   emit: versions

    script:
    """
    check_input_type.py $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

