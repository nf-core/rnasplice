process MISOPY_SETTINGS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/61/61df02ee0aef7b6fedab6906a1c054cdb2bc35dafba13d510b9c8f380708a339/data' :
        'community.wave.seqera.io/library/python_pyyaml:0610af27e7c352fd' }"

    input:
    tuple val(meta), path(bams), path(miso_data)
    val fig_width
    val fig_height

    output:
    tuple val(meta), path('miso_settings.txt'), emit: miso_settings
    path "versions.yml", topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'create_miso_settings.py'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch miso_settings.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed -n '1p' | sed 's/.*version //; s/ (.*//')"
    END_VERSIONS
    """
}
