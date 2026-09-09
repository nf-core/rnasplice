process RMATS_POST {
    tag { cond2 ? "${cond1}-${cond2}" : "${cond1}" }
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/4697fb2978d58ad7209a2d7b267bc02cf56cb8c256d787767b482c97346330d6/data'
        : 'community.wave.seqera.io/library/r-pairadise_rmats:4a32921359285283'}"

    input:
    path gtf                                     // /path/to/genome.gtf
    tuple val(contrast), val(cond1), val(meta1),  path(bam1), path(bam1_text), val(cond2), val(meta2), path(bam2), path(bam2_text), path(rmats_temp)
    val rmats_read_len                           // val params.rmats_read_len
    val rmats_splice_diff_cutoff                 // val params.rmats_splice_diff_cutoff
    val rmats_novel_splice_site                  // val params.rmats_novel_splice_site
    val rmats_min_intron_len                     // val params.rmats_min_intron_len
    val rmats_max_exon_len                       // val params.rmats_max_exon_len
    val rmats_paired_stats                       // val params.rmats_paired_stats

    output:
    path "${output_dir}/${rmats_output_dir}/*"    , emit: rmats_post
    path "${output_dir}/${logfile}"               , emit: log
    tuple val("${task.process}"), val('rmats'), eval('rmats.py --version 2>&1 | sed -e "s/v//g"'), topic: versions, emit: versions_rmats

    when:
    task.ext.when == null || task.ext.when

    script:

    // A run without a second condition still needs a named directory: Nextflow does
    // not match an output glob that starts with `./`
    output_dir = cond2 ? "${cond1}-${cond2}" : "${cond1}"
    rmats_output_dir = rmats_paired_stats ? "rmats_post_paired" : "rmats_post"

    // Only need to take meta1 as samples have same strand and read type info
    // - see rnasplice.nf input check for rmats
    def meta = meta1 instanceof List ? meta1[0] : meta1
    def args = task.ext.args ?: ''

    // Take single/paired end information from meta
    def read_type = meta.single_end ? 'single' : 'paired'

    // Default strandedness to fr-unstranded - also if user supplies "unstranded"
    def strandedness = 'fr-unstranded'
    if (meta.strandedness == 'forward') {
        strandedness  = 'fr-secondstrand'
    } else if (meta.strandedness == 'reverse') {
        strandedness  = 'fr-firststrand'
    }

    def novel_splice_sites = rmats_novel_splice_site ? '--novelSS' : ''

    // Additional args for when running with --novelSS flag
    // User defined else defauls to 50, 500
    // TODO: move these to modules.config
    def min_intron_len = ''
    def max_exon_len   = ''
    if (rmats_novel_splice_site) {
        min_intron_len = rmats_min_intron_len ? "--mil ${rmats_min_intron_len}" : '--mil 50'
        max_exon_len   = rmats_max_exon_len ? "--mel ${rmats_max_exon_len}" : '--mel 500'
    }

    def b1 = bam1_text ? "--b1 ${bam1_text}" : ''
    def b2 = ''
    def stat_flag = ''
    def paired_stats = ''

    if (cond2 && bam2_text) {
        b2 = "--b2 ${bam2_text}"
        paired_stats = rmats_paired_stats ? '--paired-stats' : ''
    }
    else {
        stat_flag = '--statoff'
    }

    logfile = rmats_paired_stats ? "rmats_post_paired.log" : "rmats_post.log"

    """
    mkdir -p ${output_dir}/${rmats_output_dir}

    rmats.py \\
        ${args} \\
        --gtf ${gtf} \\
        ${b1} \\
        ${b2} \\
        --od ${output_dir}/${rmats_output_dir} \\
        --tmp ./ \\
        -t ${read_type} \\
        --libType ${strandedness} \\
        --readLength ${rmats_read_len} \\
        --variable-read-length \\
        --nthread ${task.cpus} \\
        --tstat ${task.cpus} \\
        --cstat ${rmats_splice_diff_cutoff} \\
        --task post \\
        ${paired_stats} \\
        ${novel_splice_sites} \\
        ${min_intron_len} \\
        ${max_exon_len} \\
        ${stat_flag} \\
        --allow-clipping \\
        1> ${output_dir}/${logfile}
    """

}
