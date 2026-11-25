process CUSTOM_SPLICINGSUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/aef53fbd0e985ecd501715fa12b879f390e7d90d45f84178ee0f02d44788019c/data'
        : 'community.wave.seqera.io/library/pandas_python_pyyaml:8c93186c973f2d81'}"

    input:
    // labels[i] is a map of the kind, route and id (sample or contrast) of results[i]
    tuple val(meta), val(labels), path(results, stageAs: 'results*/*')
    path section_config // the text, headers and tool tables of the sections, assets/multiqc_splicing_sections.yml

    output:
    tuple val(meta), path("*_mqc.json"), emit: sections, optional: true
    path "versions.yml"                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'splicing_summary.py'

    stub:
    """
    touch splicing_sample_metrics_mqc.json
    touch differential_splicing_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
