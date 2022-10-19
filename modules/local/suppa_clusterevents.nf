process CLUSTEREVENTS {
    tag "$dpsi"
    label 'process_high'
    stageInMode = 'copy' 

    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path dpsi
    path psivec
    val cluster_ranges // e.g. 1-3,4-6
            
    output:
    path "*.clustvec"   , emit: clustvec 
    path "*.log"        , emit: cluster_log
    path "versions.yml" , emit: versions

            
    when:
    task.ext.when == null || task.ext.when

    script: //  Cluster events between conditions
    def st  = params.clusterevents_sigthreshold ? "-st ${params.clusterevents_sigthreshold}" : ''
    def sep = params.clusterevents_separation ? "-s ${params.clusterevents_separation}" : ''
    //def fp = file(${params.psivec})
    //def lines = fp.readLine()
    """ 
    suppa.py \\
        clusterEvents \\
        --dpsi $dpsi \\
        --psivec $psivec \\
        --dpsi-threshold ${params.clusterevents_dpsithreshold} \\
        --eps ${params.clusterevents_eps} \\
        --metric ${params.clusterevents_metric} \\
        --min-pts ${params.clusterevents_min_pts} \\
        --groups $cluster_ranges \\
        -c ${params.clusterevents_method} \\
        $st $sep -o cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
    
}