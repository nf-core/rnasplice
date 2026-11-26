process LEAFCUTTER_GROUPS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), val(samples) // samples: one list per sample, [ sample, condition, blocking variable values... ]

    output:
    tuple val(meta), path("*_groups.txt"), emit: groups

    when:
    task.ext.when == null || task.ext.when

    exec:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def groups_file = task.workDir.resolve("${prefix}_groups.txt")
    groups_file.text = samples.collect { sample -> sample.join('\t') + '\n' }.join('')
}
