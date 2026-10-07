process MERGEEVENTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/843cccb2542bb0f87d9638a058450c6f891c6af2406a42d5de9efcec3a422c18/data' :
        'community.wave.seqera.io/library/gawk:5.4.1--8dd22ef9018762ae' }"

    input:
    tuple val(meta), path(events, stageAs: 'source_ioe/*')

    output:
    tuple val(meta), path("*.ioe"), emit: ioe
    tuple val("${task.process}"), val('gawk'), eval('gawk --version | head -n 1 | sed "s/GNU Awk //; s/,.*//"'), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Keep the header of the first file only
    gawk ${args} 'FNR == 1 && NR != 1 && /^seqname/ { next } { print }' source_ioe/*.ioe > ${prefix}.ioe
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.ioe
    """
}
