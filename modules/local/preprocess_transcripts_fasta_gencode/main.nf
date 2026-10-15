process PREPROCESS_TRANSCRIPTS_FASTA_GENCODE {
    tag "${meta.id ?: fasta.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fixed.fa"), emit: fasta
    tuple val("${task.process}"), val('cut'), eval('cut --version 2>&1 | head -1 | sed "s/^.*coreutils) //"'), topic: versions, emit: versions_cut

    when:
    task.ext.when == null || task.ext.when

    script:
    def gzipped = fasta.extension == 'gz'
    def prefix = task.ext.prefix ?: (gzipped ? file(fasta.baseName).baseName : fasta.baseName)
    def command = gzipped ? 'zcat' : 'cat'
    """
    ${command} ${fasta} | cut -d "|" -f1 > ${prefix}.fixed.fa
    """

    stub:
    def gzipped = fasta.extension == 'gz'
    def prefix = task.ext.prefix ?: (gzipped ? file(fasta.baseName).baseName : fasta.baseName)
    """
    touch ${prefix}.fixed.fa
    """
}
