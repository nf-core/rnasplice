process CUSTOM_SPLICINGSUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5ad122d0a8e6da5ddadc47909d45d49a12c1232ce0d492a9592dc9e9e23e730c/data'
        : 'community.wave.seqera.io/library/pandas_python_pyyaml:599168875025982d'}"

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
