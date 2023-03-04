process CREATE_BAMLIST {
    label "process_single"

    conda "conda-forge::sed=4.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'quay.io/biocontainers/sed:4.7.0' }"

    input:
    tuple val(cond), val(meta), path(bam)

    output:
    tuple val(meta), path ("*_bamlist.txt"), emit: bam_text
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo $bam | sed 's: :,:g' > ${cond}_bamlist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -n 1 | sed 's/[^0-9.]*\\$([0-9.]*\).*/\1/')
    END_VERSIONS
    """
}
