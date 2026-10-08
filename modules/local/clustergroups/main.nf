process CLUSTERGROUPS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/459ec0a82535a33befd76060cd9d0acbeb5adef2a6eff3d07b7ee5bccb480cf1/data' :
        'community.wave.seqera.io/library/python:3.14.6--aa65261f82da0547' }"

    input:
    tuple val(meta), path(psivec)

    output:
    tuple val(meta), path("*_groups.txt"), emit: groups
    tuple val(meta), eval('cat *_groups.txt'), emit: ranges
    tuple val("${task.process}"), val('python'), eval('python3 -c "import platform; print(platform.python_version())"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - <<'EOF'
    import itertools

    with open("${psivec}") as handle:
        samples = handle.readline().rstrip("\\n").split("\\t")

    # Trim the replicate number, e.g. GBR_1 -> GBR
    conditions = [sample.rsplit("_", 1)[0] for sample in samples]

    ranges = []
    start = 1
    for _condition, group in itertools.groupby(conditions):
        end = start + len(list(group)) - 1
        ranges.append(f"{start}-{end}")
        start = end + 1

    if len(ranges) != 2:
        raise ValueError("Column numbers have to be continuous, with no overlapping or missing columns between them.")

    with open("${prefix}_groups.txt", "w") as handle:
        handle.write(",".join(ranges) + "\\n")
    EOF
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_groups.txt
    """
}
