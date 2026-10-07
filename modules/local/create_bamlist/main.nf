process CREATE_BAMLIST {
    tag "${contrast}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2d2d3fbc9bc617edaba31d8c03a305397730a6ad06070318c242938e399428b0/data' :
        'community.wave.seqera.io/library/sed:4.10--711bc92627aef086' }"

    input:
    tuple val(contrast), val(cond1), path(bam1), val(cond2), path(bam2)

    output:
    tuple val(contrast), path("${cond1}_bamlist.txt"), emit: bam_list1
    tuple val(contrast), path("${cond2}_bamlist.txt"), emit: bam_list2, optional: true
    tuple val("${task.process}"), val('sed'), eval('sed --version 2>&1 | head -n1 | sed "s/.* //"'), topic: versions, emit: versions_sed

    when:
    task.ext.when == null || task.ext.when

    script:
    // rMATS reads the BAM files of a condition as one comma separated line. A single
    // condition run has no second condition, so its bam list 2 is not written
    def b1_list = bam1 instanceof List ? bam1.join(' ') : bam1
    def b2_list = bam2 instanceof List ? bam2.join(' ') : bam2
    def bam2_cmd = (cond2 && b2_list) ? "echo ${b2_list} | sed 's: :,:g' > ${cond2}_bamlist.txt" : ''
    """
    echo ${b1_list} | sed 's: :,:g' > ${cond1}_bamlist.txt
    ${bam2_cmd}
    """

    stub:
    def b2_list = bam2 instanceof List ? bam2.join(' ') : bam2
    def bam2_cmd = (cond2 && b2_list) ? "touch ${cond2}_bamlist.txt" : ''
    """
    touch ${cond1}_bamlist.txt
    ${bam2_cmd}
    """
}
