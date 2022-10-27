process EDGER_EXON {
    tag "$samplesheet"
    label "process_single"

    conda (params.enable_conda ? "bioconda::bioconductor-edger=3.36.0 conda-forge::r-statmod=1.4.36" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0' :
        'quay.io/biocontainers/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0' }"

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
    run_edger_exon.R featurecounts $samplesheet $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-edger:  \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """

}
