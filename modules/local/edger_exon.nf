process EDGER_EXON {
    tag "$samplesheet"
    label "process_low"

    conda (params.enable_conda ? "bioconda::bioconductor-edger=3.36.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-edger:3.36.0--r41hc247a5b_2' :
        'quay.io/biocontainers/bioconductor-edger:3.36.0--r41hc247a5b_2' }"

    input:
    path ("featurecounts/*")
    path samplesheet

    output:
    path "DGEList.rds"  , emit: edger_exon_dge
    path "DGEGLM.rds"   , emit: edger_exon_glm
    path "DGELRT.*.rds" , emit: edger_exon_lrt
    path "*.csv"        , emit: edger_exon_csv
    path "versions.yml" , emit: versions

    script:
    """
    run_edger_exon.R featurecounts $samplesheet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-edger:  \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """

}