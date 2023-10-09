process DEXSEQ_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::htseq=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py310ha14a713_0' :
        'biocontainers/htseq:2.0.2--py310ha14a713_0' }"

    input:
    tuple val(meta), path(bam), path (gff)
    val alignment_quality                   // val params.alignment_quality

    output:
    tuple val(meta), path("*.clean.count.txt"), emit: dexseq_clean_txt
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_type = meta.single_end ? '' : '-p yes'

    def alignment_quality = "-a ${alignment_quality}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '-s yes'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '-s reverse'
    } else if (meta.strandedness == 'unstranded') {
        strandedness = '-s no'
    }

    """
    dexseq_count.py $gff $read_type -f bam $bam -r pos ${prefix}.clean.count.txt $alignment_quality $strandedness

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('htseq').version)")
    END_VERSIONS
    """
}
