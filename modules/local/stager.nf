process STAGER {
    label "process_medium"

    conda "bioconda::bioconductor-stager=1.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-stager:1.12.0--r40_0' :
        'quay.io/biocontainers/bioconductor-stager:1.12.0--r40_0' }"

    input:
    tuple path(feature_rds), path(gene_rds)    // dxr.rds or drimseq_d.rds || qvals.rds or res.txp.rds
    val analysis_type    // dexseq or drimseq

    output:
    path "stageRTx.*.rds"              , emit: stager_rds
    path "getAdjustedPValues.*.rds"    , emit: stager_padj_rds
    path "getAdjustedPValues.*.tsv"    , emit: stager_padj_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    run_stager.R $feature_rds $gene_rds $analysis_type

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-stager: \$(Rscript -e "library(stageR); cat(as.character(packageVersion('stageR')))")
    END_VERSIONS
    """
}
