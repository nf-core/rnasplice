process CREATE_BAMLIST_SINGLE {
    label "process_single"

    conda "conda-forge::sed=4.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'quay.io/biocontainers/sed:4.7.0' }"

    input:
    tuple val(cond), val(meta), path(bam)

    output:
    tuple val(cond), val(meta), path(bam), path("${cond}_bamlist.txt"), emit: bam_text
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo $bam | sed 's: :,:g' > ${cond}_bamlist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
