process DIFFSPLICE {
    tag "${cond1}-${cond2}"
    label 'process_high'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path events
    tuple val(cond1), val(cond2), path(tpm1), path(tpm2)
    tuple val(cond1), val(cond2), path(psi1), path(psi2)
    val prefix

    output:
    tuple val(cond1), val(cond2), path("*.dpsi")       , emit: dpsi
    tuple val(cond1), val(cond2), path("*.psivec")     , emit: psivec
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: //  Calculate differential analysis between conditions

    def diffsplice_method = params.diffsplice_method ?: 'empirical' // empirical or classical

    def gc     = params.diffsplice_gene_correction ? "-gc" : ''
    def pa     = params.diffsplice_paired ? "-pa" : ''
    def median = params.diffsplice_median ? "-me" : ''

    def diffsplice_area          = params.diffsplice_area ?: '1000'       // default 1000
    def diffsplice_lower_bound   = params.diffsplice_lower_bound ?: '0'   // default 0
    def diffsplice_alpha         = params.diffsplice_alpha ?: '0.05'      // default 0.05
    def diffsplice_tpm_threshold = params.diffsplice_tpm_threshold ?: '0' // default 0
    def diffsplice_nan_threshold = params.diffsplice_nan_threshold ?: '0' // default 0

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
