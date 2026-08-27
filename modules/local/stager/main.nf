process STAGER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/41a73be6620eaa0bb67220f2b3e71168e4724c09556c76b65c2ad88f383d7bd4/data':
        'community.wave.seqera.io/library/bioconductor-stager:1.32.0--854b4a08a81d122d' }"

    input:
    tuple val(meta), path(feature_tsv), path(gene_tsv)
    val analysis_type

    output:
    tuple val(meta), path("stageRTx.*.rds")          , emit: stager_rds
    tuple val(meta), path("getAdjustedPValues.*.rds"), emit: stager_padj_rds
    tuple val(meta), path("getAdjustedPValues.*.tsv"), emit: stager_padj_tsv
    path "versions.yml"                              , topic: versions, emit: versions_stager

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'run_stager.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    echo ${args}

    touch stageRTx.${prefix}.rds
    touch getAdjustedPValues.${prefix}.rds
    touch getAdjustedPValues.${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stageR: \$(Rscript -e "library(stageR); cat(as.character(packageVersion('stageR')))")
    END_VERSIONS
    """
}
