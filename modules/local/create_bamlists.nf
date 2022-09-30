process CREATE_BAMLISTS {
    //tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"
    
    input:
    tuple val(cond), val(meta), path(bam1)
    tuple val(cond), val(meta), path(bam2)

    output:
    path("bamlist_group1.txt"), emit:bam_group1
    path("bamlist_group2.txt"), emit:bam_group2

    exec:
    task.workDir.resolve('bamlist_group1.txt').text = bam1.join(',')
    task.workDir.resolve('bamlist_group2.txt').text = bam2.join(',')
}
