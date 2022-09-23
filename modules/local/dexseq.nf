process DEXSEQ {
    label "process_medium"

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0 bioconda::bioconductor-drimseq=1.18.0 bioconda::bioconductor-stager=1.12.0" : null)
    container "docker.io/yuukiiwa/nanoseq:dexseq"
    // need a multitool container for r-base, dexseq, stager, drimseq and on quay hub

    input:
    path drimseq_filter_rds

    output:
    path "dxd.rds"                , emit: dexseq_rds
    path "dxr.rds"                , emit: dexseq_results_rds
    path "dxr.tsv"                , emit: dexseq_results_tsv
    path "qval.rds"               , emit: qval_rds
    path "dxr.g.tsv"              , emit: dexseq_results_q_tsv
    path "versions.yml"           , emit: versions

    script:
    """
    run_dexseq.R $drimseq_filter_rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}