process LEAFCUTTER_GTFTOEXONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/leafcutter:2.0.3--pyhd8ed1ab_0'
        : 'quay.io/biocontainers/leafcutter:2.0.3--pyhd8ed1ab_0'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.exons.txt.gz"), emit: exons
    // leafcutter has no option to print its version, so it is read from the Python package
    tuple val("${task.process}"), val('leafcutter'), eval("python -c \"import importlib.metadata as m; print(m.version('leafcutter'))\""), topic: versions, emit: versions_leafcutter

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    leafcutter-gtf-to-exons \\
        ${gtf} \\
        ${prefix}.exons.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.exons.txt.gz
    """
}
