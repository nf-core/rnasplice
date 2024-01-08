process CLUSTEREVENTS {
    tag "${cond1}-${cond2}"
    label 'process_high'
    stageInMode = 'copy'

    conda "bioconda::suppa=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa:2.3--py36_0' :
        'biocontainers/suppa:2.3--py36_0' }"

    input:
    tuple val(cond1), val(cond2), path(dpsi)
    tuple val(cond1), val(cond2), path(psivec)
    val group_ranges // e.g. 1-3,4-6
    val prefix
    val clusterevents_dpsithreshold    // val params.clusterevents_dpsithreshold
    val clusterevents_eps              // val params.clusterevents_eps
    val clusterevents_metric           // val params.clusterevents_metric
    val clusterevents_min_pts          // val params.clusterevents_min_pts
    val clusterevents_method           // val params.clusterevents_method
    val clusterevents_sigthreshold     // val params.clusterevents_sigthreshold
    val clusterevents_separation       // val params.clusterevents_separation

    output:
    path "*.clustvec"   , emit: clustvec
    path "*.log"        , emit: cluster_log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: //  Cluster events between conditions

    def clusterevents_sigthreshold  = clusterevents_sigthreshold ? "-st ${params.clusterevents_sigthreshold}" : ''
    def clusterevents_separation = clusterevents_separation ? "-s ${params.clusterevents_separation}" : ''

    """
    suppa.py \\
        clusterEvents \\
        --dpsi $dpsi \\
        --psivec $psivec \\
        --dpsi-threshold $clusterevents_dpsithreshold \\
        --eps $clusterevents_eps \\
        --metric $clusterevents_metric \\
        --min-pts $clusterevents_min_pts \\
        --groups $group_ranges \\
        --clustering $clusterevents_method \\
        $clusterevents_sigthreshold $clusterevents_separation -o ${cond1}-${cond2}_${prefix}_cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """

}
