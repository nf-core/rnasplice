process SUBREAD_FLATTENGTF {
    tag "$annotation"
    label 'process_single'

    conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'biocontainers/subread:2.0.1--hed695b0_0' }"

    input:
    path(annotation)

    output:
    path "annotation.saf", emit: saf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    flattenGTF $args -a $annotation -o annotation.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(flattenGTF -v 2>&1) | sed -e "s/flattenGTF v//g")
    END_VERSIONS
    """
}
