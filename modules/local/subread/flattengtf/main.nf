process SUBREAD_FLATTENGTF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0'
        : 'biocontainers/subread:2.0.1--hed695b0_0'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.saf"), emit: saf
    tuple val("${task.process}"), val('subread'), eval('flattenGTF -v 2>&1 | head -1 | sed -e "s/flattenGTF v//g"'), topic: versions, emit: versions_subread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flattenGTF \\
        $args \\
        -a $gtf \\
        -o ${prefix}.saf
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}.saf
    """
}
