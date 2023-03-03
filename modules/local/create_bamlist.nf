process CREATE_BAMLIST {
    //tag "$meta.id"
    label "process_low"

    conda "conda-forge::sed=4.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'quay.io/biocontainers/sed:4.7.0' }"

    input:
    tuple val(cond), val(meta), path(bam)

    output:
    path("*_bamlist.txt"), emit: bam_text

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo $bam | sed 's: :,:g' > ${cond}_bamlist.txt
    """
}
