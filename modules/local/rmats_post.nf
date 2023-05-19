process RMATS_POST {
    tag "$cond1-$cond2"
    label 'process_high'

    conda 'bioconda::r-pairadise=1.0.0 bioconda::rmats=4.1.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae:121e48ab4817ec619c157a346458efca1ccf3c0a-0' :
        'quay.io/biocontainers/mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae:121e48ab4817ec619c157a346458efca1ccf3c0a-0' }"

    input:
    path gtf                                     // /path/to/genome.gtf
    path bam_group1                              // path("bamlist_group1.txt")
    path bam_group2                              // path("bamlist_group2.txt")
    tuple val(cond1), val(meta1), path(bam1)     // [condition1, [condition1_metas], [condition1_bams]]
    tuple val(cond2), val(meta2), path(bam2)     // [condition2, [condition2_metas], [condition2_bams]]
    path ("$cond1-$cond2/rmats_temp/*")          // rmats temp folder from rmats prep step

    output:
    path "$cond1-$cond2/rmats_post/*"        , emit: rmats_post
    path "$cond1-$cond2/rmats_post.log"      , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Only need to take meta1 as samples have same strand and read type info
    // - see rnasplice.nf input check for rmats
    def meta = meta1[0]
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "$cond1-$cond2"

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

    // Take read length input as user defined else defaults to 40
    def read_len = params.rmats_read_len ?: '40'

    // User defined label if samples are paired and paired stats required
    // Rmats uses bioconductor-pairadise Rscript downstream
    def paired_stats = params.rmats_paired_stats ? '--paired-stats' : ''

    // Take read length input as user defined else default to 0.0001
    def splice_cutoff = params.rmats_splice_diff_cutoff ?: '0.0001'

    // Whether user wants to run with novel splice sites flag
    def novel_splice_sites = params.rmats_novel_splice_site ? '--novelSS' : ''

    // Additional args for when running with --novelSS flag
    // User defined else defauls to 50, 500
    def min_intron_len = ''
    def max_exon_len   = ''
    if (params.rmats_novel_splice_site) {
        min_intron_len = params.rmats_min_intron_len ? "--mil ${params.rmats_min_intron_len}" : '--mil 50'
        max_exon_len   = params.rmats_max_exon_len ? "--mel ${params.rmats_max_exon_len}" : '--mel 500'
    }

    """
    mkdir -p $prefix/rmats_post

    rmats.py \\
        --gtf $gtf \\
        --b1 $bam_group1 \\
        --b2 $bam_group2 \\
        --od $prefix/rmats_post \\
        --tmp $prefix/rmats_temp \\
        -t $read_type \\
        --libType $strandedness \\
        --readLength $read_len \\
        --variable-read-length \\
        --nthread $task.cpus \\
        --tstat $task.cpus \\
        --cstat $splice_cutoff \\
        --task post \\
        $paired_stats \\
        $novel_splice_sites \\
        $min_intron_len \\
        $max_exon_len \\
        --allow-clipping \\
        1> $prefix/rmats_post.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(echo \$(rmats.py --version) | sed -e "s/v//g")
    END_VERSIONS
    """

}
