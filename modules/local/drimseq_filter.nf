process DRIMSEQ_FILTER {
    label "process_medium"

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0 bioconda::bioconductor-drimseq=1.18.0 bioconda::bioconductor-stager=1.12.0" : null)
    container "docker.io/yuukiiwa/nanoseq:dexseq"
    // need a multitool container for r-base, dexseq, stager, drimseq and on quay hub

    input:
    path txi         // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    path tx2gene     // path: *.tx2gene.tsv
    path samplesheet // path: /path/to/samplesheet.csv

    output:
    path "*.rds"                , emit: drimseq_filter_rds
    path "*.csv"                , emit: drimseq_filter_csv
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: '' 
    """
    run_drimseq.r $txi $tx2gene $samplesheet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-drimseq: \$(Rscript -e "library(DRIMSeq); cat(as.character(packageVersion('DRIMSeq')))")
    END_VERSIONS
    """
}