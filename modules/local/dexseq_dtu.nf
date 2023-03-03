process DEXSEQ_DTU {
    label "process_high"

    conda "bioconda::bioconductor-dexseq=1.36.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.36.0--r40_0' :
        'quay.io/biocontainers/bioconductor-dexseq:1.36.0--r40_0' }"

    input:
    path drimseq_sample_data
    path drimseq_d_counts

    output:
    path "dxd.rds"                , emit: dexseq_rds
    path "dxr.rds"                , emit: dexseq_results_rds
    path "dxr.tsv"                , emit: dexseq_results_tsv
    path "qval.rds"               , emit: qval_rds
    path "dxr.g.tsv"              , emit: dexseq_results_q_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def denominator = params.dtu_lfc_denominator ?: ""

    """
    run_dexseq_dtu.R $drimseq_sample_data $drimseq_d_counts ${task.cpus} $denominator

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}
