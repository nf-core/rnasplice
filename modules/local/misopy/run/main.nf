process MISOPY_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a6a585c9b85aee50b8da33d7d4b27d209174ea7c7279a2ac73ceda6b8d4d193/data':
        'community.wave.seqera.io/library/python_misopy:9d5a390611c447f5' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(miso_index)

    output:
    tuple val(meta), path("miso_data/*")    , emit: miso
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('misopy'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\'misopy\').version)"'), topic: versions, emit: versions_misopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    miso --run \\
        ${miso_index} \\
        ${bam} \\
        --output-dir miso_data/${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p miso_data/${prefix}
    touch miso_data/${prefix}/miso_results
    """
}
