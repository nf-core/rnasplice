process CREATE_BAMLIST {
    label 'process_single'

    conda "conda-forge::sed=4.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'biocontainers/sed:4.7.0' }"

    input:
    tuple val(cond1), val(meta1), path(bam1)
    tuple val(cond2), val(meta2), path(bam2)

    output:
    tuple val(cond1), val(meta1), path(bam1), path("${cond1}_bamlist.txt"), emit: bam1_text
    tuple val(cond2), val(meta2), path(bam2), path("${cond2}_bamlist.txt"), emit: bam2_text
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo $bam1 | sed 's: :,:g' > ${cond1}_bamlist.txt
    echo $bam2 | sed 's: :,:g' > ${cond2}_bamlist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
