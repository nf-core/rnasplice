process SPLIT_FILES {
    tag "$tpm_psi"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0':
        'quay.io/biocontainers/r-base:3.4.2' }"

    input:
    path tpm_psi
    path samplesheet
    val output_type  // either .tpm or .psi
    val calc_ranges  // true/false calculate ranges
    val prefix       // output file prefix - transcript or local

    output:
    path "*.tpm"        , optional : true , emit: tpms
    path "*.psi"        , optional : true , emit: psis
    path "ranges.txt"   , optional : true , emit: ranges
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    suppa_split_file.R \\
        $tpm_psi \\
        $samplesheet \\
        $output_type \\
        $calc_ranges \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
