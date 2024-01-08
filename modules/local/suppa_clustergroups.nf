process CLUSTERGROUPS {
    tag "${cond1}-${cond2}"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(cond1), val(cond2), path(psivec)

    output:
    stdout

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    suppa_groups.py $psivec
    """

}
