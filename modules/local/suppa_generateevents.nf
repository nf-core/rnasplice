process GENERATE_EVENTS {
    tag "$gtf"
    label 'process_low'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'biocontainers/suppa:2.3--py36_0' }"

    input:
    path gtf
    val file_type
    val generateevents_boundary       // val params.generateevents_boundary
    val generateevents_threshold      // val params.generateevents_threshold
    val generateevents_exon_length    // val params.generateevents_exon_length
    val generateevents_event_type     // val params.generateevents_event_type
    val generateevents_pool_genes     // val params.generateevents_pool_genes

    output:
    path "events.*"     , emit: events
    path "versions.yml" , emit: versions
    path "events_*.*"   , emit: eventstype, optional : true
    // Declaring as optional as these are produced only in local events and not transcript events

    when:
    task.ext.when == null || task.ext.when

    script:

    // If pool_genes is set to true then include the -p parameter
    def poolgenes = generateevents_pool_genes ? "-p" : ''

    // Calculate local events and combine all the ioe events files
    if (file_type == 'ioe') {
        """
        suppa.py \\
            generateEvents \\
            -i $gtf \\
            -f $file_type \\
            -o events \\
            -e $generateevents_event_type \\
            -b $generateevents_boundary \\
            -t $generateevents_threshold \\
            -l $generateevents_exon_length \\
            $poolgenes

        awk 'FNR==1 && NR!=1 { while (/^seqname/) getline; }  1 {print}' *.ioe > events.ioe

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            suppa: \$(pip show suppa| sed -e '/Version/!d'| sed 's/Version: //g')
        END_VERSIONS
        """
    }
    // Calculate transcript events
    else if (file_type == 'ioi') {
        """
        suppa.py \\
            generateEvents \\
            -i $gtf \\
            -f $file_type \\
            -o events \\
            $poolgenes

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
        END_VERSIONS
        """
    }
}
