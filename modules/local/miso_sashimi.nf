process MISO_SASHIMI {
    label 'process_single'

    conda "conda-forge::python=2.7 bioconda::misopy=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27h9801fc8_5' :
        'biocontainers/misopy:0.5.4--py27h9801fc8_5' }"

    input:
    path index_path
    tuple path(miso_settings), val(miso_gene)
    path ("bam_files/*")   // Need bam and bai in working directory
    path ("miso_data/*")   // Need Miso files in working directory

    output:
    path "sashimi/*"           , emit: sashimi
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sashimi_plot --plot-event $miso_gene $index_path $miso_settings --output-dir sashimi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed "s/Python //g")
        misopy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('misopy').version)")
    END_VERSIONS
    """

}
