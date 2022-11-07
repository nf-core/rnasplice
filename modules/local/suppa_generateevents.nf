process GENERATE_EVENTS {
    tag "$gtf"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path gtf
    val file_type

    output:
    path "events.*"     , emit: events
    path "versions.yml" , emit: versions
    path "events_*.*"   , emit: eventstype, optional : true
    /* Declaring as optional as these are produced only in local events and not transcript events */

    when:
    task.ext.when == null || task.ext.when

    script:

    // Calculate the AS events using GTF

    // Incase of Local Events get the list of events from nextflow.config
    def list_events = params.local_events ?: 'SE SS MX RI FL' // default all possible

    // If pool_genes is set to true then include the -p parameter
    def poolgenes = params.pool_genes ? "-p" : ''

    // Calculate local events and combine all the ioe events files
    if (file_type == 'ioe') {
        """
        suppa.py \\
            generateEvents \\
            -i $gtf \\
            -f $file_type \\
            -o events \\
            -e $list_events \\
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
