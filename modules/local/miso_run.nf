process MISO_RUN {
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::python=2.7 bioconda::misopy=0.5.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27h9801fc8_5' :
        'quay.io/biocontainers/misopy:0.5.4--py27h9801fc8_5' }"

    input:
    path miso_index
    tuple val(meta), path(bams), path(bais)
    val miso_read_len

    output:
    path "miso_data/*"         , emit: miso
    path "versions.yml"        , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    miso --run ${miso_index} $bams --output-dir miso_data/${prefix} --read-len $miso_read_len

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        misopy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('misopy').version)")
    END_VERSIONS
    """

}
