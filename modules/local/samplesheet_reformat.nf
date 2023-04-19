process SAMPLESHEET_REFORMAT {
    tag "$params.input"
    label 'process_single'

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

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

