process SPLIT_FILES {
    tag "${tpm_psi}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7aa202fe46bc3ec5eb02d21c659ee8d93b1995724712fface3837047691651a2/data':
        'community.wave.seqera.io/library/r-base:4.6.1--e4f1a108384f0df8' }"

    input:
    tuple val(meta), path(tpm_psi)
    path samplesheet
    val output_type  // either .tpm or .psi
    val calc_ranges  // true/false calculate ranges

    output:
    path "*.tpm"        , optional : true , emit: tpms
    path "*.psi"        , optional : true , emit: psis
    path "*_ranges.txt" , optional : true , emit: ranges
    path "versions.yml", topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'suppa_split_file.R'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    echo ${args}

    touch ${prefix}.tpm
    touch ${prefix}.psi
    touch ${prefix}_ranges.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa_split_file: "\$(Rscript --version 2>&1 | sed -n '1p' | sed 's/.*version //; s/ (.*//')"
    END_VERSIONS
    """
}
