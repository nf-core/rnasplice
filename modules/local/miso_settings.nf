process MISO_SETTINGS {
    label 'process_single'

    conda "conda-forge::parsimonious"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/parsimonious:0.10.0' :
        'biocontainers/parsimonious:0.10.0' }"

    input:
    path miso
    path bams
    val fig_width
    val fig_height

    output:
    path 'miso_settings.txt'          , emit: miso_settings
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    create_miso_settings.py \\
        $args \\
        --bams $bams \\
        --name $miso \\
        --width $fig_width \\
        --height $fig_height \\
        --output 'miso_settings.txt'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsimonious: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('parsimonious').version)")
    END_VERSIONS
    """

}
