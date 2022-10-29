process DRIMSEQ_FILTER {
    label "process_medium"

    conda (params.enable_conda ? "bioconda::bioconductor-drimseq=1.18.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-drimseq:1.18.0--r40_0' :
        'quay.io/biocontainers/bioconductor-drimseq:1.18.0--r40_0' }"

    input:
    path txi                  // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    path tximport_tx2gene     // path: tximport.tx2gene.tsv
    path samplesheet          // path: /path/to/samplesheet.csv

    output:
    path "d.rds"                , emit: drimseq_filter_rds
    path "sample.data.tsv"      , emit: drimseq_sample_data
    path "d.counts.tsv"         , emit: drimseq_d_counts
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''

    def min_samps_gene_expr = params.min_samps_gene_expr ?: 0
    def min_samps_feature_expr = params.min_samps_feature_expr ?: 0
    def min_samps_feature_prop = params.min_samps_feature_prop ?: 0
    def min_feature_expr = params.min_feature_expr ?: 0
    def min_feature_prop = params.min_feature_prop ?: 0
    def min_gene_expr = params.min_gene_expr ?: 0

    """
    run_drimseq_filter.R $txi $tximport_tx2gene $samplesheet \\
        $min_samps_gene_expr \\
        $min_samps_feature_expr \\
        $min_samps_feature_prop \\
        $min_feature_expr \\
        $min_feature_prop \\
        $min_gene_expr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-drimseq: \$(Rscript -e "library(DRIMSeq); cat(as.character(packageVersion('DRIMSeq')))")
    END_VERSIONS
    """
}
