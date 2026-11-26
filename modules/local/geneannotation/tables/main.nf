process GENEANNOTATION_TABLES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/459ec0a82535a33befd76060cd9d0acbeb5adef2a6eff3d07b7ee5bccb480cf1/data' :
        'community.wave.seqera.io/library/python:3.14.6--aa65261f82da0547' }"

    input:
    tuple val(meta), path(tx2gene) // tx2gene: `gffread --table` of transcript_id, gene_id, the extra attributes (the first being the gene name), gene_biotype and gene_type

    output:
    tuple val(meta), path("*.gene_annotation.tsv")      , emit: gene
    tuple val(meta), path("*.transcript_annotation.tsv"), emit: transcript
    path "versions.yml"                                , topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'geneannotation_tables.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'gene_id\\tgene_name\\tgene_biotype\\n' > ${prefix}.gene_annotation.tsv
    printf 'transcript_id\\tgene_id\\tgene_name\\tgene_biotype\\n' > ${prefix}.transcript_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //g')"
    END_VERSIONS
    """
}
