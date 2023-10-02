process PSIPERISOFORM {
    tag "$tpm"
    label 'process_medium'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'biocontainers/suppa:2.3--py36_0' }"

    input:
    path gtf
    path tpm

    output:
    path "suppa_isoform.psi"  , emit: psi
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // Calculate the psi values isoform level
    """
    suppa.py \\
        psiPerIsoform \\
        -g $gtf \\
        -e $tpm \\
        -o suppa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}
