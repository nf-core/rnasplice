process EDGER_EXON {
    label "process_high"

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0 bioconda::bioconductor-drimseq=1.18.0 bioconda::bioconductor-stager=1.12.0" : null)
    container "quay.io/biocontainers/bioconductor-edger"
    // need a multitool container for r-base, dexseq, stager, drimseq and on quay hub

    input:
    path ("dexseq_clean_counts/*")     // path dexseq_clean_counts
    path samplesheet                   // path: samplesheet

    output:
    path "dxd_exon.rds"           , emit: dexseq_exon_rds
    path "versions.yml"           , emit: versions

    script:
    """
    run_edger_exon.R {} {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}