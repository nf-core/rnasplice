process DEXSEQ_EXON {
    label "process_high"

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0 bioconda::bioconductor-drimseq=1.18.0 bioconda::bioconductor-stager=1.12.0" : null)
    container "docker.io/yuukiiwa/nanoseq:dexseq"
    // need a multitool container for r-base, dexseq, stager, drimseq and on quay hub

    input:
    path ("dexseq_clean_counts/*")     // path dexseq_clean_counts
    path gff                           // path dexseq_gff
    path samplesheet                   // path samplesheet
    val read_method                    // htseq or featurecounts

    output:
    path "dxd_exon.rds"           , emit: dexseq_exon_rds
    path "dxr_exon.rds"           , emit: dexseq_exon_results_rds
    path "dxr_exon.tsv"           , emit: dexseq_exon_results_tsv
    path "qval_exon.rds"          , emit: qval_exon_rds
    path "dxr_exon.g.tsv"         , emit: dexseq_exon_results_q_tsv
    path "versions.yml"           , emit: versions

    script:
    def demoninator = params.lfc_denominator ?: ""

    """
    run_dexseq_exon.R dexseq_clean_counts $gff $samplesheet $read_method ${task.cpus} $demoninator

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}