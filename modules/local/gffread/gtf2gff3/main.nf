process GFFREAD_GTF2GFF3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4':
        'biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3
    tuple val("${task.process}"), val('gffread'), eval("gffread --version"), topic: versions, emit: versions_gffread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-L --keep-genes'
    def prefix = task.ext.prefix ?: "${meta.id}"
    // `transcript` features are renamed to `mRNA` so that the resulting GFF3 is
    // accepted by downstream tools (e.g. MISO) that follow the GFF3 spec strictly.
    """
    gffread \\
        $gtf \\
        $args \\
        | awk -F'\\t' -vOFS='\\t' '{ gsub("transcript", "mRNA", \$3); print}' \\
        > ${prefix}.gff3
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    """
}
