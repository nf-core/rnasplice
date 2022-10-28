process DIFFSPLICE {
    tag "$tpms"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py_2' }"

    input:
    path events
    path tpms
    path psis
    val prefix

    output:
    path "*.dpsi"       , emit: dpsi
    path "*.psivec"     , emit: psivec
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
        -p $psis \\
        -e $tpms \\
        -o ${prefix}_diffsplice

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
}
