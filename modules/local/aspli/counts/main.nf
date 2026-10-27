process ASPLI_COUNTS {
    tag "${meta.id}"
    // Only binGenome can use more than one core, and gbCounts and jCounts run on one,
    // so the task asks for a single cpu. ASpli reads every BAM file into memory in
    // full, so memory grows with the size of the BAM files and is set in modules.config.
    // These resources have not been measured on a real dataset yet
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-aspli:2.16.0--r44hdfd78af_0'
        : 'biocontainers/bioconductor-aspli:2.16.0--r44hdfd78af_0'}"

    input:
    tuple val(meta), path(targets), path(bams), path(bais) // targets: TSV with sample, condition, bam, single_end and strandedness columns
    tuple val(meta2), path(gtf)
    val read_length                                       // integer: --minReadLength of ASpli, 0 to estimate it from the BAM files

    output:
    tuple val(meta), path("*.features.rds")      , emit: features
    tuple val(meta), path("*.gbcounts.rds")      , emit: gbcounts
    tuple val(meta), path("*.jcounts.rds")       , emit: jcounts
    tuple val(meta), path("*_counts.tsv")        , emit: counts
    tuple val(meta), path("*.{pir,pju,psi}.tsv") , emit: junction_usage
    tuple val(meta), path("*.settings.tsv")      , emit: settings
    path "versions.yml"                          , topic: versions, emit: versions_aspli

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'aspli_counts.R'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    for object in features gbcounts jcounts; do
        touch ${prefix}.\${object}.rds
    done
    for table in gene bin junction e1i ie2; do
        touch ${prefix}.\${table}_counts.tsv
    done
    for table in junctions.pir junctions.pju ir.pir alt.psi es.psi; do
        touch ${prefix}.\${table}.tsv
    done
    touch ${prefix}.settings.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*\$//')
        bioconductor-aspli: \$(Rscript -e "cat(as.character(packageVersion('ASpli')))")
    END_VERSIONS
    """
}
