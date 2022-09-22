process DEXSEQ_ANNOTATION {
    tag "$gtf"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0' }"
	
    input:
    path gtf

    output:
    path "*.gff"        , emit:gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    dexseq_prepare_annotation.py $gtf genome.gff
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    	r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq:  \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
	python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}