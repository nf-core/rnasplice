process CLUSTEREVENTS {
    tag "$dpsi"
    label 'process_high'
    stageInMode = 'copy'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path dpsi
    path psivec
    val cluster_ranges // e.g. 1-3,4-6
    val prefix

    output:
    path "*.clustvec"   , emit: clustvec
    path "*.log"        , emit: cluster_log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: //  Cluster events between conditions
    def st  = params.clusterevents_sigthreshold ? "-st ${params.clusterevents_sigthreshold}" : ''
    def sep = params.clusterevents_separation ? "-s ${params.clusterevents_separation}" : ''

    def clusterevents_dpsithreshold = params.clusterevents_dpsithreshold ?: '0.05'  // default 0.05
    def clusterevents_eps           = params.clusterevents_eps ?: '0.05'            // default 0.05
    def clusterevents_metric        = params.clusterevents_metric ?: 'euclidean'    // default euclidean
    def clusterevents_min_pts       = params.clusterevents_min_pts ?: '20'          // default 20
    def clusterevents_method        = params.clusterevents_method ?: 'DBSCAN'       // default DBSCAN

    """

    suppa.py \\
        clusterEvents \\
        --dpsi $dpsi \\
        --psivec $psivec \\
        --dpsi-threshold $clusterevents_dpsithreshold \\
        --eps $clusterevents_eps \\
        --metric $clusterevents_metric \\
        --min-pts $clusterevents_min_pts \\
        --groups $cluster_ranges \\
        --clustering $clusterevents_method \\
        $st $sep -o ${prefix}_cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """

}
