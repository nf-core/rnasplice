process DRIMSEQ_FILTER {
    label 'process_medium'

    conda "bioconda::bioconductor-drimseq=1.18.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-drimseq:1.18.0--r40_0' :
        'biocontainers/bioconductor-drimseq:1.18.0--r40_0' }"

    input:
    path txi                    // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    path tximport_tx2gene       // path: tximport.tx2gene.tsv
    path samplesheet            // path: /path/to/samplesheet.csv
    val min_samps_gene_expr     // val params.min_samps_gene_expr
    val min_samps_feature_expr  // val params.min_samps_feature_expr
    val min_samps_feature_prop  // val params.min_samps_feature_prop
    val min_feature_expr        // val params.min_feature_expr
    val min_feature_prop        // val params.min_feature_prop
    val min_gene_expr           // val params.min_gene_expr

    output:
    path "dmDSdata.rds"  , emit: drimseq_dataset_rds
    path "samples.tsv"   , emit: drimseq_samples_tsv
    path "counts.tsv"    , emit: drimseq_counts_tsv
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

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
