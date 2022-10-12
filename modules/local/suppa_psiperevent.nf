process PSIPEREVENT {

    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path ioi_ioe
    path tpm
    val type

    output:
    path "*.psi"        , emit: psi
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // Calculate the psi values of events
    """
    suppa.py \\
        psiPerEvent \\
        -i $ioi_ioe \\
        -e $tpm \\
        -o $type
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}