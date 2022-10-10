process DIFFSPLICE {
  
    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py_2' }"

    input:
    path events
    path tpms
    path psis
        
    output:
    path "versions.yml", emit: versions
    path "*.dpsi" , emit: dpsi 
    path "*.psivec" , emit: psivec

    when:
    task.ext.when == null || task.ext.when

    script: //  Calculate differential analysis between conditions
    def gc = params.diffsplice_gene_correction ? "-gc" : ''
    def pa = params.diffsplice_paired ? "-pa" : ''
    def median = params.diffsplice_median ? "-me" : ''
        
    """
    suppa.py \\
        diffSplice \\
        -m ${params.diffsplice_method} \\
        $gc $pa -s -c $median \\
        -a ${params.diffsplice_area} \\
        -l ${params.diffsplice_lower_bound} \\
        -al ${params.diffsplice_alpha} \\
        -th ${params.diffsplice_tpm_threshold} \\
        -nan ${params.diffsplice_nan_threshold} \\
        -i $events \\
        -p $psis \\
        -e $tpms \\
        -o diffsplice
     
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}