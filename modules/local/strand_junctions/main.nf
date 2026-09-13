process STRAND_JUNCTIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52399c97123d96b3b6eca039a907a3536c9c6aaddcf5d09c8ae50ec82df35e95/data' :
        'community.wave.seqera.io/library/python_pyyaml:0610af27e7c352fd' }"

    input:
    tuple val(meta), path(junc)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("stranded/*.junc"), emit: junc
    path "versions.yml"                     , topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'strand_junctions.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir stranded
    touch stranded/${prefix}.junc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //g')"
    END_VERSIONS
    """
}
