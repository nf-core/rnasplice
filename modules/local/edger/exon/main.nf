process EDGER_EXON {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0'
        : 'biocontainers/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0'}"

    input:
    tuple val(meta) , path(feature_counts, stageAs: "featurecounts/*") // path: <sample>.featureCounts.tsv
    tuple val(meta2), path(samplesheet)                               // path: /path/to/samplesheet.csv
    tuple val(meta3), path(contrastsheet)                             // path: /path/to/contrastsheet.csv

    output:
    tuple val(meta), path("DGEList.rds")  , emit: edger_exon_dge
    tuple val(meta), path("DGEGLM.rds")   , emit: edger_exon_glm
    tuple val(meta), path("DGELRT.*.rds") , emit: edger_exon_lrt
    tuple val(meta), path("*.csv")        , emit: edger_exon_csv
    tuple val(meta), path("*.pdf")        , emit: edger_exon_pdf
    path "versions.yml"                   , topic: versions, emit: versions_edger

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'run_edger_exon.R'

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}

    touch DGEList.rds
    touch DGEGLM.rds
    touch DGELRT.exprs.rds
    touch DGELRT.usage.rds

    for contrast in \$(tail -n +2 ${contrastsheet} | cut -d ',' -f 1); do
        touch "contrast_\${contrast}.exprs.csv"
        for test in simes gene exon; do
            touch "contrast_\${contrast}.usage.\${test}.csv"
            touch "contrast_\${contrast}.usage.\${test}.pdf"
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*\$//')
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}
