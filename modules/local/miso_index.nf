process MISO_INDEX {
    label 'process_high'

    conda "conda-forge::python=2.7 bioconda::misopy=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27h9801fc8_5' :
        'biocontainers/misopy:0.5.4--py27h9801fc8_5' }"

    input:
    path gff3
    val index

    output:
    path index           , emit: miso_index
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    index_gff --index $gff3 $index
    parse_miso_index.py -p $index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed "s/Python //g")
        misopy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('misopy').version)")
    END_VERSIONS
    """
}
