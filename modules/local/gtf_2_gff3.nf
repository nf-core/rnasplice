process GTF_2_GFF3 {

    label "process_single"
    conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'quay.io/biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    path gtf

    output:
    path "*.gff3"           , emit: gff3
    path "versions.yml"     , emit: versions

    script:
    """
    gffread $gtf -L | awk -F'\t' -vOFS='\t' '{ gsub("transcript", "mRNA", \$3); print}' > ${gtf.baseName}_genes.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
