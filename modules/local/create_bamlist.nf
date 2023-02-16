process CREATE_BAMLIST {
    //tag "$meta.id"
    label "process_low"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(cond), val(meta), path(bam)

    output:
    path("*_bamlist.txt"), emit:bam_text

    script:
    """
    echo $bam | sed 's: :,:g' > ${cond}_bamlist.txt
    """
}
