process CREATEBAMLIST {
    tag "${contrast}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65d91e147d41e8367a773f5587cf71b53a18af4eb8494feb8a4e18f423184be3/data' :
        'community.wave.seqera.io/library/sed:4.9--da997413f41d23b0' }"

    input:
    tuple val(contrast), val(cond1), val(ids1), path(bam1), val(cond2), val(ids2), path(bam2)

    output:
    tuple val(contrast), path("${cond1}_bamlist.txt"), emit: bam_list1
    tuple val(contrast), path("${cond2}_bamlist.txt"), emit: bam_list2, optional: true
    tuple val(contrast), path("${task.ext.prefix ?: contrast}.samples.tsv"), emit: samples
    tuple val("${task.process}"), val('sed'), eval('sed --version 2>&1 | head -n1 | sed "s/.* //"'), topic: versions, emit: versions_sed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: contrast
    // rMATS reads the BAM files of a condition as one comma separated line. A single
    // condition run has no second condition, so its bam list 2 is not written
    def bams1 = bam1 instanceof List ? bam1 : [bam1]
    def bams2 = bam2 instanceof List ? bam2 : [bam2]
    def bam2_cmd = (cond2 && bams2) ? "echo ${bams2.join(' ')} | sed 's: :,:g' > ${cond2}_bamlist.txt" : ''
    // rMATS never names its samples: the n-th value of a `_SAMPLE_1` or `_SAMPLE_2` column
    // is the n-th BAM file of the first or second bam list, so the same order is tabulated
    def rows1 = [ids1, bams1].transpose().withIndex().collect { pair, idx ->
        [contrast, 'SAMPLE_1', idx + 1, pair[0], cond1, cond2 ? 'target' : 'profiling', pair[1]].join('\\t')
    }
    def rows2 = !cond2 ? [] : [ids2, bams2].transpose().withIndex().collect { pair, idx ->
        [contrast, 'SAMPLE_2', idx + 1, pair[0], cond2, 'reference', pair[1]].join('\\t')
    }
    def sample_rows = rows1 + rows2
    """
    echo ${bams1.join(' ')} | sed 's: :,:g' > ${cond1}_bamlist.txt
    ${bam2_cmd}

    printf 'contrast\\trmats_group\\tposition\\tsample\\tcondition\\trole\\tbam\\n' > ${prefix}.samples.tsv
    ${sample_rows.collect { row -> "printf '${row}\\n' >> ${prefix}.samples.tsv" }.join('\n    ')}
    """

    stub:
    def prefix = task.ext.prefix ?: contrast
    def b2_list = bam2 instanceof List ? bam2.join(' ') : bam2
    def bam2_cmd = (cond2 && b2_list) ? "touch ${cond2}_bamlist.txt" : ''
    """
    touch ${cond1}_bamlist.txt
    ${bam2_cmd}
    touch ${prefix}.samples.tsv
    """
}
