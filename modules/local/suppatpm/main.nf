process SUPPATPM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(tpm)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: tpm
    tuple val("${task.process}"), val('cut'), eval('cut --version 2>&1 | head -1 | sed "s/^.*coreutils) //"'), topic: versions, emit: versions_cut

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // SUPPA reads the first column as transcript IDs only when the header has one field
    // fewer than the rows, so the ID column name is dropped along with the gene column
    """
    head -n 1 ${tpm} | cut -f 3- > ${prefix}.txt
    tail -n +2 ${tpm} | cut -f 1,3- >> ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
