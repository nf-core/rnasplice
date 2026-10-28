process ASPLI_DU {
    tag "${meta.id}"
    // These resources have not been measured on a real dataset yet
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-aspli:2.16.0--r44hdfd78af_0'
        : 'biocontainers/bioconductor-aspli:2.16.0--r44hdfd78af_0'}"

    input:
    tuple val(meta), path(features), path(gbcounts), path(jcounts) // meta: [ id: contrast, treatment, control ]

    output:
    tuple val(meta), path("*.{genes_de,bins_du}.tsv")                             , emit: gb_du
    tuple val(meta), path("*.{locale,anchor,ir,es,alt}_{junctions,clusters}.tsv") , emit: j_du
    tuple val(meta), path("*.splicing_report_*.tsv")                              , emit: splicing_report
    tuple val(meta), path("*.integrated_signals.tsv")                             , emit: integrated_signals
    tuple val(meta), path("*_integrated_signals")                                 , emit: integrated_signals_html
    tuple val(meta), path("*.aspli_du.rds")                                       , emit: rds
    path "versions.yml"                                                           , topic: versions, emit: versions_aspli

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'aspli_du.R'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    for table in genes_de bins_du locale_junctions locale_clusters anchor_junctions anchor_clusters ir_junctions es_junctions alt_junctions splicing_report_bins splicing_report_locale splicing_report_anchor integrated_signals; do
        touch ${prefix}.\${table}.tsv
    done
    mkdir ${prefix}_integrated_signals
    touch ${prefix}_integrated_signals/integratedSignals.html
    touch ${prefix}.aspli_du.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*\$//')
        bioconductor-aspli: \$(Rscript -e "cat(as.character(packageVersion('ASpli')))")
    END_VERSIONS
    """
}
