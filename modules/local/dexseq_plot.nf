process DEXSEQ_PLOT {
    label "process_single"

    conda (params.enable_conda ? "bioconda::bioconductor-dexseq=1.36.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.36.0--r40_0' :
        'quay.io/biocontainers/bioconductor-dexseq:1.36.0--r40_0' }"

    input:
    path dxr                           // path dxr_exon.rds
    val n_dexseq_plot                  // number to plot

    output:
    path "dexseq_plot.pdf"  , emit: dexseq_plot_pdf
    path "versions.yml"     , emit: versions

    script:
    """
    run_dexseq_plot.R $dxr $n_dexseq_plot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}
