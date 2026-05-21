process MULTIQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.18--pyhdfd78af_0' :
        'biocontainers/multiqc:1.18--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), emit: plots, optional: true
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .
    """

    stub:
    """
    touch multiqc_data
    touch multiqc_plots
    touch multiqc_report.html
    """
}
