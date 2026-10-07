process DEXSEQ_EXON {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/673dc7e5ec7cbb244b3464de0b78558cec773fb24184c60980aaae6bbfeadbcf/data' :
        'community.wave.seqera.io/library/bioconductor-dexseq:1.56.0--52df10363b3f1635' }"

    input:
    path ("dexseq_clean_counts/*")    // path: dexseq_clean_counts
    path gff                          // path: dexseq_gff
    path samplesheet                  // path: samplesheet
    path contrastsheet                // path: contrastsheet
    val ntop                          // val: n_dexseq_plot

    output:
    path "DEXSeqDataSet.*.rds"  , emit: dexseq_exon_dataset_rds
    path "DEXSeqResults.*.rds"  , emit: dexseq_exon_results_rds
    path "perGeneQValue.*.rds"  , emit: dexseq_gene_results_rds
    path "DEXSeqResults.*.csv"  , emit: dexseq_exon_results_csv
    path "perGeneQValue.*.csv"  , emit: dexseq_gene_results_csv
    path "plotDEXSeq.*.pdf"     , emit: dexseq_plot_results_pdf
    tuple val("${task.process}"), val('r-base'), eval('R --version 2>&1 | head -n 1 | sed "s/^.*version //; s/ .*$//"'), topic: versions, emit: versions_R
    tuple val("${task.process}"), val('bioconductor-dexseq'), eval('Rscript -e "library(DEXSeq); cat(as.character(packageVersion(\'DEXSeq\')))"'), topic: versions, emit: versions_dexseq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_dexseq_exon.R dexseq_clean_counts $gff $samplesheet $contrastsheet $ntop ${task.cpus}
    """

    stub:
    """
    # `run_dexseq_exon.R` names its output files after the `contrast` column of the
    # contrastsheet, which it finds by its header
    for contrast in \$(awk -F ',' 'NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "contrast") col = i; next } { print \$col }' ${contrastsheet}); do
        touch DEXSeqDataSet.\${contrast}.rds
        touch DEXSeqResults.\${contrast}.rds
        touch perGeneQValue.\${contrast}.rds
        touch DEXSeqResults.\${contrast}.csv
        touch perGeneQValue.\${contrast}.csv
        touch plotDEXSeq.\${contrast}.pdf
    done
    """
}
