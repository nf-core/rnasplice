process DEXSEQ_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::htseq=2.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    	'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py310ha14a713_0' :
	'quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0' }"
	
    input:
    path gff
    tuple val(meta), path(bam)

    output:
    path "*.txt"        , emit: dexseq_txt
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read_type = meta.single_end ? '' : '-p yes'    
    def alignment_quality = "-a ${params.alignment_quality}"
    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '-s yes'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '-s reverse'
    } else if (meta.strandedness == 'unstranded') {
	strandedness = '-s no'
    }

    """
    dexseq_count.py $gff $read_type -f bam $bam -r pos count.txt $alignment_quality $strandedness
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	htseq: \$(pip show htseq | sed -e '/Version/!d'| sed 's/Version: //g')
    END_VERSIONS
    """
}