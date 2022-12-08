process DEXSEQ_EXON {
    label "process_high"

    conda (params.enable_conda ? "bioconda::bioconductor-dexseq=1.36.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dexseq:1.36.0--r40_0' :
        'quay.io/biocontainers/bioconductor-dexseq:1.36.0--r40_0' }"

    input:
    path ("dexseq_clean_counts/*")     // path dexseq_clean_counts
    path gff                           // path dexseq_gff
    path gtf                           // path featurecounts_gtf
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
    def denominator = params.deu_lfc_denominator ?: ""

    """
    run_dexseq_exon.R dexseq_clean_counts $gff $gtf $samplesheet $read_method ${task.cpus} $denominator

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}
