process MISO_RUN {
    label 'process_medium'

    conda "conda-forge::python=2.7 bioconda::misopy=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27h9801fc8_5' :
        'biocontainers/misopy:0.5.4--py27h9801fc8_5' }"

    input:
    path miso_index
    tuple val(meta), path(bams), path(bais)
    val miso_read_len

    output:
    tuple val(meta), path("miso_data/*")    , emit: miso
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    miso --run ${miso_index} $bams --output-dir miso_data/${meta.id} --read-len $miso_read_len

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        misopy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('misopy').version)")
    END_VERSIONS
    """

}
