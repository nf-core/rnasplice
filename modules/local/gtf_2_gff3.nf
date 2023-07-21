process GTF_2_GFF3 {
    label 'process_single'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_3' :
        'biocontainers/gffread:0.12.7--hdcf5f25_3' }"

    input:
    path gtf

    output:
    path "*.gff3"           , emit: gff3
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread $gtf -L --keep-genes | awk -F'\\t' -vOFS='\\t' '{ gsub("transcript", "mRNA", \$3); print}' > ${gtf.baseName}_genes.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
