process CLUSTERGROUPS {
    tag "${cond1}-${cond2}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52399c97123d96b3b6eca039a907a3536c9c6aaddcf5d09c8ae50ec82df35e95/data' :
        'community.wave.seqera.io/library/python_pyyaml:0610af27e7c352fd' }"

    input:
    tuple val(meta), val(cond1), val(cond2), path(psivec)

    output:
    tuple val(meta), val(cond1), val(cond2), path("*_groups.txt"), emit: groups
    path "versions.yml", topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'cluster_groups.py'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${cond1}-${cond2}"
    """
    echo ${args}

    touch ${prefix}_groups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa_groups: "\$(python3 --version 2>&1 | sed -n '1p' | sed 's/.*version //; s/ (.*//')"
    END_VERSIONS
    """
}
