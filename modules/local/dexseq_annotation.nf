process DEXSEQ_ANNOTATION {
    tag "$gtf"
    label 'process_medium'

    conda "bioconda::htseq=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py310ha14a713_0' :
    'biocontainers/htseq:2.0.2--py310ha14a713_0' }"

    input:
    path gtf         // path gtf file
    val aggregation  // val params.aggregation

    output:
    path "*.gff"        , emit: gff
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "DEXSeq"

    def aggregation = aggregation ? '' : '-r no'

    """
    dexseq_prepare_annotation.py $gtf ${prefix}.gff $aggregation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('htseq').version)")
    END_VERSIONS
    """
}
