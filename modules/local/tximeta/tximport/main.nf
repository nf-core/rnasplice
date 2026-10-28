process TXIMETA_TXIMPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.28.2--r45hdfd78af_0':
        'biocontainers/bioconductor-tximeta:1.28.2--r45hdfd78af_0' }"

    input:
    tuple val(meta), path("salmon/*")
    path tx2gene

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

    path "versions.yml"                         , topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix            = task.ext.prefix ?: meta.id
    def ignore_tx_version = task.ext.args ?: 'false'
    template "tximport.R"

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.txi.rds
    touch ${prefix}.txi.s.rds
    touch ${prefix}.txi.ls.rds
    touch ${prefix}.txi.dtu.rds

    touch ${prefix}.gi.rds
    touch ${prefix}.gi.s.rds
    touch ${prefix}.gi.ls.rds

    touch ${prefix}.gene_tpm.tsv
    touch ${prefix}.gene_counts.tsv
    touch ${prefix}.gene_tpm_scaled.tsv
    touch ${prefix}.gene_counts_scaled.tsv
    touch ${prefix}.gene_tpm_length_scaled.tsv
    touch ${prefix}.gene_counts_length_scaled.tsv

    touch ${prefix}.transcript_tpm.tsv
    touch ${prefix}.transcript_counts.tsv
    touch ${prefix}.transcript_tpm_scaled.tsv
    touch ${prefix}.transcript_counts_scaled.tsv
    touch ${prefix}.transcript_tpm_length_scaled.tsv
    touch ${prefix}.transcript_counts_length_scaled.tsv
    touch ${prefix}.transcript_tpm_dtu_scaled.tsv
    touch ${prefix}.transcript_counts_dtu_scaled.tsv

    touch tximport.tx2gene.tsv

    touch suppa_tpm.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e 'cat(as.character(packageVersion("tximeta")))')
        bioconductor-tximport: \$(Rscript -e 'cat(as.character(packageVersion("tximport")))')
    END_VERSIONS
    """
}
