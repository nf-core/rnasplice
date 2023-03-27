process MISO_SASHIMI {
    label "process_single"

    conda (params.enable_conda ? "conda-forge::python=2.7 bioconda::misopy=0.5.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27h9801fc8_5' :
        'quay.io/biocontainers/misopy:0.5.4--py27h9801fc8_5' }"

    input:
    path miso_path
    tuple path(miso_settings), val(miso_gene)

    output:
    path "sashimi/*"           , emit: sashimi
    path "versions.yml"        , emit: versions


    script:
    """
    sashimi_plot --plot-event $miso_gene $miso_path $miso_settings --output-dir sashimi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed "s/Python //g")
        misopy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('misopy').version)")
    END_VERSIONS
    """

}
