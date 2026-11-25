process CUSTOM_SPLICINGSUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6bb45b3a3fed2ebd55af1322a927ea11a751b8b62e542021dbf5be0d375abbb7/data'
        : 'community.wave.seqera.io/library/python:3.12.12--2411dfb0536aecdf'}"

    input:
    // labels[i] is a map of the kind, route and id (sample or contrast) of results[i]
    tuple val(meta), val(labels), path(results, stageAs: 'results*/*')

    output:
    tuple val(meta), path("*_mqc.json"), emit: sections, optional: true
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'splicing_summary.py'

    stub:
    """
    touch splicing_sample_metrics_mqc.json
    touch differential_splicing_mqc.json
    """
}
