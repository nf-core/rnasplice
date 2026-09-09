process RMATS_PREP {
    tag { cond2 ? "${cond1}-${cond2}" : "${cond1}" }
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8b/8b7a2184a4c9e054a13811c086f2d6095637be0f519191db74c35650085f5c20/data' :
        'community.wave.seqera.io/library/rmats:4.3.0--177f3a2035a879e5' }"

    input:
    path gtf                                     // /path/to/genome.gtf
    tuple val(contrast),  val(cond1), val(meta1), path(bam1), path(bam1_text), val(cond2),  val(meta2), path(bam2), path(bam2_text)
    val rmats_read_len                           // val params.rmats_read_len
    val rmats_splice_diff_cutoff                 // val params.rmats_splice_diff_cutoff
    val rmats_novel_splice_site                  // val params.rmats_novel_splice_site
    val rmats_min_intron_len                     // val params.rmats_min_intron_len
    val rmats_max_exon_len                       // val params.rmats_max_exon_len

    output:
    tuple val(contrast), path("${output_dir}/rmats_temp/*") , emit: rmats_temp
    path "${output_dir}/rmats_prep.log"                     , emit: log
    tuple val("${task.process}"), val('rmats'), eval('rmats.py --version 2>&1 | sed -e "s/v//g"'), topic: versions, emit: versions_rmats

    when:
    task.ext.when == null || task.ext.when

    script:
    // A run without a second condition still needs a named directory: Nextflow does
    // not match an output glob that starts with `./`
    output_dir = cond2 ? "${cond1}-${cond2}" : "${cond1}"

     // Only need to take meta1 as samples have same strand and read type info
    // Only need to take meta1 as samples have same strand and read type info
    // - see rnasplice.nf input check for rmats
    def meta = meta1 instanceof List ? meta1[0] : meta1
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: output_dir

    // Take single/paired end information from meta
    def read_type = meta.single_end ? 'single' : 'paired'

    // Default strandedness to fr-unstranded - also if user supplies "unstranded"
    def strandedness = 'fr-unstranded'

    // Change strandedness based on user samplesheet input
    if (meta.strandedness == 'forward') {
        strandedness  = 'fr-secondstrand'
    } else if (meta.strandedness == 'reverse') {
        strandedness  = 'fr-firststrand'
    }

    // Whether user wants to run with novel splice sites flag
    def novel_splice_sites = rmats_novel_splice_site ? '--novelSS' : ''
    def min_intron_len = ''
    def max_exon_len   = ''
    if (rmats_novel_splice_site) {
        min_intron_len = rmats_min_intron_len ? "--mil ${rmats_min_intron_len}" : '--mil 50'
        max_exon_len   = rmats_max_exon_len ? "--mel ${rmats_max_exon_len}" : '--mel 500'
    }

    def b1 = bam1_text ? "--b1 ${bam1_text}" : ''
    def b2 = (cond2 && bam2_text) ? "--b2 ${bam2_text}" : ''

    """
    mkdir -p ${prefix}/rmats_temp
    mkdir -p ${prefix}/rmats_prep

    rmats.py \\
        ${args} \\
        ${b1} \\
        ${b2} \\
        -t ${read_type} \\
        --libType ${strandedness} \\
        --nthread ${task.cpus} \\
        --tstat ${task.cpus} \\
        --gtf ${gtf} \\
        --allow-clipping \\
        --readLength ${rmats_read_len} \\
        --variable-read-length \\
        --cstat ${rmats_splice_diff_cutoff} \\
        --task prep \\
        ${novel_splice_sites} \\
        ${min_intron_len} \\
        ${max_exon_len} \\
        --od ${prefix}/rmats_prep \\
        --tmp ${prefix}/rmats_temp 1> ${prefix}/rmats_prep.log
    """

}
