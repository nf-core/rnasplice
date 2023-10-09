process TXIMPORT {
    tag "$tx2gene"
    label 'process_medium'

    conda "bioconda::bioconductor-tximeta=1.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0' :
        'biocontainers/bioconductor-tximeta:1.8.0--r40_0' }"

    input:

    path ("salmon/*")
    path tx2gene // path to tx2gene.tsv

    output:

    path "*txi.rds"                              , emit: txi
    path "*txi.s.rds"                            , emit: txi_s
    path "*txi.ls.rds"                           , emit: txi_ls
    path "*txi.dtu.rds"                          , emit: txi_dtu

    path "*gi.rds"                               , emit: gi
    path "*gi.s.rds"                             , emit: gi_s
    path "*gi.ls.rds"                            , emit: gi_ls

    path "*gene_tpm.tsv"                        , emit: tpm_gene
    path "*gene_counts.tsv"                     , emit: counts_gene
    path "*gene_tpm_scaled.tsv"                 , emit: tpm_gene_scaled
    path "*gene_counts_scaled.tsv"              , emit: counts_gene_scaled
    path "*gene_tpm_length_scaled.tsv"          , emit: tpm_gene_length_scaled
    path "*gene_counts_length_scaled.tsv"       , emit: counts_gene_length_scaled

    path "*transcript_tpm.tsv"                  , emit: tpm_transcript
    path "*transcript_counts.tsv"               , emit: counts_transcript
    path "*transcript_tpm_scaled.tsv"           , emit: tpm_transcript_scaled
    path "*transcript_counts_scaled.tsv"        , emit: counts_transcript_scaled
    path "*transcript_tpm_length_scaled.tsv"    , emit: tpm_transcript_length_scaled
    path "*transcript_counts_length_scaled.tsv" , emit: counts_transcript_length_scaled
    path "*transcript_tpm_dtu_scaled.tsv"       , emit: tpm_transcript_dtu_scaled
    path "*transcript_counts_dtu_scaled.tsv"    , emit: counts_transcript_dtu_scaled

    path "tximport.tx2gene.tsv"                 , emit: tximport_tx2gene

    path "suppa_tpm.txt"                        , emit: suppa_tpm

    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnasplice/bin/
    """
    tximport.R $tx2gene salmon salmon.merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
