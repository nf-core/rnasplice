process RMATS_POST {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/4697fb2978d58ad7209a2d7b267bc02cf56cb8c256d787767b482c97346330d6/data'
        : 'community.wave.seqera.io/library/r-pairadise_rmats:4a32921359285283'}"

    input:
    tuple val(meta), path(rmats, stageAs: 'rmats_tmp/*'), path(bam_list1), path(bam_list2)
    tuple val(meta2), path(gtf)
    val read_length

    output:
    tuple val(meta), path("${prefix}/*.MATS.{JC,JCEC}.txt")     , emit: mats
    tuple val(meta), path("${prefix}/fromGTF.*.txt")            , emit: from_gtf
    tuple val(meta), path("${prefix}/{JC,JCEC}.raw.input.*.txt"), emit: raw_input
    tuple val(meta), path("${prefix}/summary.txt")              , emit: summary
    tuple val(meta), path("${prefix}.log")                      , emit: log
    tuple val("${task.process}"), val('rmats'), eval('rmats.py --version | sed -e "s/v//g"'), topic: versions, emit: versions_rmats
    tuple val("${task.process}"), val('pairadise'), eval('Rscript -e "cat(as.character(packageVersion(\'PAIRADISE\')))"'), topic: versions, emit: versions_pairadise

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // A single condition has nothing to test against, so rMATS only counts the events
    def b2   = bam_list2 ? "--b2 ${bam_list2}" : ''
    def stat = bam_list2 ? '' : '--statoff'
    """
    mkdir -p ${prefix}

    rmats.py \\
        --task post \\
        ${args} \\
        --gtf ${gtf} \\
        --b1 ${bam_list1} \\
        ${b2} \\
        ${stat} \\
        --readLength ${read_length} \\
        --nthread ${task.cpus} \\
        --tstat ${task.cpus} \\
        --tmp rmats_tmp \\
        --od ${prefix} \\
        > ${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}
    for event in SE MXE A3SS A5SS RI; do
        touch ${prefix}/\${event}.MATS.JC.txt
        touch ${prefix}/\${event}.MATS.JCEC.txt
        touch ${prefix}/fromGTF.\${event}.txt
        touch ${prefix}/fromGTF.novelJunction.\${event}.txt
        touch ${prefix}/fromGTF.novelSpliceSite.\${event}.txt
        touch ${prefix}/JC.raw.input.\${event}.txt
        touch ${prefix}/JCEC.raw.input.\${event}.txt
    done
    touch ${prefix}/summary.txt
    touch ${prefix}.log
    """
}
