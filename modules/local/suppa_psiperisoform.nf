process PSIPERISOFORM {
  
    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path gtf
    path tpm
 
    output:
    path "*.psi"  , emit: psi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // Calculate the psi values
    """
    suppa.py \\
        psiPerIsoform \\
        -g $gtf \\
        -e $tpm \\
        -o gtf_tpm
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}