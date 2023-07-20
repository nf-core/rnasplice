process PSIPEREVENT {
    tag "$tpm"
    label 'process_medium'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'biocontainers/suppa:2.3--py36_0' }"

    input:
    path ioe
    path tpm
    val psiperevent_total_filter   // val params.psiperevent_total_filter

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
        -f $psiperevent_total_filter \\
        -o suppa_local

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}
