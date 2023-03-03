process RMATS_PREP_SINGLE {
    label "process_medium"

    conda 'bioconda::r-pairadise=1.0.0 bioconda::rmats=4.1.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c5a18b683684b2b5cb608213a267834373186743:1c831c11c46ebe0c54bcb67918bc35fedf1c43ee-0' :
        'quay.io/biocontainers/mulled-v2-c5a18b683684b2b5cb608213a267834373186743:1c831c11c46ebe0c54bcb67918bc35fedf1c43ee-0' }"

    input:
    path gtf                                     // /path/to/genome.gtf
    path bam_group1                              // path("bamlist_group1.txt")
    tuple val(cond1), val(meta1), path(bams)     // [condition1, [condition1_metas], [condition1_bams]]

    output:
    path "rmats_temp/*"      , emit: rmats_temp
    path "rmats_prep.log"    , emit: log
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Only need to take meta1 as samples have same strand and read type info
    // - see rnasplice.nf input check for rmats
    def meta = meta1[0]
    def args = task.ext.args ?: ''

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
    rmats.py \\
        --b1 $bam_group1 \\
        -t $read_type \\
        --libType $strandedness \\
        --nthread $task.cpus \\
        --gtf $gtf \\
        --allow-clipping \\
        --readLength $read_len \\
        --variable-read-length \\
        --cstat $splice_cutoff \\
        --task prep \\
        $novel_splice_sites \\
        $min_intron_len \\
        $max_exon_len \\
        --tmp rmats_temp \\
        --od rmats_prep 1> rmats_prep.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(echo \$(rmats.py --version) | sed -e "s/v//g")
    END_VERSIONS
    """
    }
