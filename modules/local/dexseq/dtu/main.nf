process DEXSEQ_DTU {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/673dc7e5ec7cbb244b3464de0b78558cec773fb24184c60980aaae6bbfeadbcf/data' :
        'community.wave.seqera.io/library/bioconductor-dexseq:1.56.0--52df10363b3f1635' }"

    input:
    path drimseq_sample_data
    path drimseq_d_counts
    path drimseq_contrast_data

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "DEXSeqResults.*.tsv"  , emit: dexseq_exon_results_tsv
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "perGeneQValue.*.tsv"  , emit: dexseq_gene_results_tsv
    tuple val("${task.process}"), val('r-base'), eval('R --version 2>&1 | head -n 1 | sed "s/^.*version //; s/ .*$//"'), topic: versions, emit: versions_R
    tuple val("${task.process}"), val('bioconductor-dexseq'), eval('Rscript -e "library(DEXSeq); cat(as.character(packageVersion(\'DEXSeq\')))"'), topic: versions, emit: versions_dexseq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_dexseq_dtu.R ${drimseq_sample_data} \\
        ${drimseq_contrast_data} \\
        ${drimseq_d_counts}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}

    # `run_dexseq_dtu.R` names its output files after the contrast, which it builds by
    # joining the treatment and the control column of the contrast table with a dash
    for contrast in \$(tail -n +2 ${drimseq_contrast_data} | awk -F ',' '{ print \$2 "-" \$3 }'); do
        touch DEXSeqDataSet.\${contrast}.rds
        touch DEXSeqResults.\${contrast}.rds
        touch DEXSeqResults.\${contrast}.tsv
        touch perGeneQValue.\${contrast}.rds
        touch perGeneQValue.\${contrast}.tsv
    done
    """
}
