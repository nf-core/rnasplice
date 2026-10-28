process DRIMSEQ_DMFILTER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-drimseq:1.18.0--r40_0':
        'biocontainers/bioconductor-drimseq:1.18.0--r40_0' }"

    input:
    path txi               // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    path tximport_tx2gene  // path: tximport.tx2gene.tsv
    path samplesheet       // path: /path/to/samplesheet.csv

    output:
    path "dmDSdata.rds", emit: drimseq_dataset_rds
    path "samples.tsv" , emit: drimseq_samples_tsv
    path "counts.tsv"  , emit: drimseq_counts_tsv
    path "versions.yml", topic: versions, emit: versions_drimseq

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'run_drimseq_filter.R'

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}

    touch dmDSdata.rds
    touch samples.tsv
    touch counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*\$//')
        bioconductor-drimseq: \$(Rscript -e "library(DRIMSeq); cat(as.character(packageVersion('DRIMSeq')))")
    END_VERSIONS
    """
}
