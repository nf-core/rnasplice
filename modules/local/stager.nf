process STAGER {
    label "process_medium"

    conda "bioconda::bioconductor-stager=1.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-stager:1.12.0--r40_0' :
        'quay.io/biocontainers/bioconductor-stager:1.12.0--r40_0' }"

    input:

    path rds          // dxr.rds or drimseq_d.rds
    val analysis_type // dexseq or drimseq
    path results_rds  // qvals.rds or res.txp.rds

    output:

    path "*.stageR.padj.tsv"    , emit: stager_padj_tsv
    path "*.stageR.padj.rds"    , emit: stager_padj_rds
    path "*.stageRObj.rds"      , emit: stager_rds
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    run_stager.R $rds $analysis_type $results_rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-stager: \$(Rscript -e "library(stageR); cat(as.character(packageVersion('stageR')))")
    END_VERSIONS
    """
}
