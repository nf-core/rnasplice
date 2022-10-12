process CLUSTEREVENTS {
  
    conda (params.enable_conda ? "bioconda::suppa" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/suppa%3A2.3--py36_0' :
        'quay.io/biocontainers/suppa:2.3--py36_0' }"

    input:
    path dpsi
    path psivec
            
    output:
    path "versions.yml", emit: versions
    path "*.clustvec" , emit: clustvec 
    path "*.log" , emit: cluster_log
 
            
    when:
    task.ext.when == null || task.ext.when

    script: //  Cluster events between conditions
    """
    cp $dpsi dpsi
    
     suppa.py \\
        clusterEvents \\
        --dpsi dpsi \\
        --psivec $psivec \\
        --sig-threshold 0.05 \\
        --eps 0.2 \\
        --separation 0.11 \\
        -dt 0.2 --min-pts 10 \\
        --groups 1-3,4-6 \\
        -c OPTICS -o cluster
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('suppa').version)")
    END_VERSIONS
    """
    
}