process PSIPEREVENT {
    tag "$tpm"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path ioe
    path tpm

    output:
    path "suppa_local.psi"    , emit: psi
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // Calculate the psi values of local events
    """
    suppa.py \\
        psiPerEvent \\
        -i $ioe \\
        -e $tpm \\
        -o suppa_local

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}