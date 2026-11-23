process GENEANNOTATION_ANNOTATE {
    tag "${table.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/459ec0a82535a33befd76060cd9d0acbeb5adef2a6eff3d07b7ee5bccb480cf1/data' :
        'community.wave.seqera.io/library/python:3.14.6--aa65261f82da0547' }"

    input:
    tuple val(meta), path(table)
    tuple val(meta2), path(gene_annotation)
    val id_column

    output:
    tuple val(meta), path("*.annotated.{tsv,csv}"), emit: annotated
    path "versions.yml"                           , topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'geneannotation_annotate.py'

    stub:
    // The annotated copy of a CSV table is a CSV, and of any other table a TSV
    def prefix = task.ext.prefix ?: table.baseName
    def extension = table.extension == 'csv' ? 'csv' : 'tsv'
    """
    touch ${prefix}.annotated.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //g')"
    END_VERSIONS
    """
}
