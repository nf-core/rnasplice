/**
 * Reads the read lengths that rMATS ran with from the scripts of its tasks, as rMATS
 * does not record them in its outputs.
 */
class RmatsReadLength {

    /**
     * Maps each rMATS step ('prep' and 'post') that ran under the work directory to the
     * distinct `--readLength` values of its tasks. Stub tasks run no rMATS command and
     * are not counted.
     */
    static Map<String, List<String>> fromWorkDir(Object work_dir) {
        def read_lengths = [:]
        new File(work_dir.toString()).eachFileRecurse { script ->
            if (script.name != '.command.sh') {
                return
            }
            def text = script.text
            def step = text =~ /rmats\.py\s+\\?\s*--task (\S+)/
            def read_len = text =~ /--readLength (\S+)/
            if (step.find() && read_len.find()) {
                read_lengths.computeIfAbsent(step.group(1)) { [] as Set } << read_len.group(1)
            }
        }
        read_lengths.collectEntries { step_name, values -> [ step_name, values.sort() ] }
    }

}
