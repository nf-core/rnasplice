process DIFFSPLICE {
    tag "${cond1}-${cond2}"
    label 'process_high'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'biocontainers/suppa:2.3--py36_0' }"

    input:
    path events
    tuple val(cond1), val(cond2), path(tpm1), path(tpm2)
    tuple val(cond1), val(cond2), path(psi1), path(psi2)
    val prefix
    val diffsplice_method             // val params.diffsplice_method
    val diffsplice_area               // val params.diffsplice_area
    val diffsplice_lower_bound        // val params.diffsplice_lower_bound
    val diffsplice_alpha              // val params.diffsplice_alpha
    val diffsplice_tpm_threshold      // val params.diffsplice_tpm_threshold
    val diffsplice_nan_threshold      // val params.diffsplice_nan_threshold
    val diffsplice_gene_correction    // val params.diffsplice_gene_correction
    val diffsplice_paired             // val params.diffsplice_paired
    val diffsplice_median             // val params.diffsplice_median

    output:
    tuple val(cond1), val(cond2), path("*.dpsi")       , emit: dpsi
    tuple val(cond1), val(cond2), path("*.psivec")     , emit: psivec
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: //  Calculate differential analysis between conditions

    def gc     = diffsplice_gene_correction ? "-gc" : ''
    def pa     = diffsplice_paired ? "-pa" : ''
    def median = diffsplice_median ? "-me" : ''

    """
    suppa.py \\
        diffSplice \\
        -m $diffsplice_method \\
        $gc $pa -s -c $median \\
        -a $diffsplice_area \\
        -l $diffsplice_lower_bound \\
        -al $diffsplice_alpha \\
        -th $diffsplice_tpm_threshold \\
        -nan $diffsplice_nan_threshold \\
        -i $events \\
        -p $psi1 $psi2 \\
        -e $tpm1 $tpm2 \\
        -o ${cond1}-${cond2}_${prefix}_diffsplice

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}
