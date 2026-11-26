/**
 * Reads the tables leafcutter-ds writes for snapshots that do not depend on the machine.
 *
 * The statistics come from a numerical fit whose result differs from one machine to
 * another past the second or third significant digit, so they cannot be snapshotted even
 * rounded. Only the columns that identify each cluster or intron and the outcome of its
 * test are kept: the cluster or intron id, the test status, the degrees of freedom and
 * the genes.
 */
class LeafcutterDsTable {

    static final List<String> STABLE_COLUMNS = ['cluster', 'intron', 'status', 'df', 'genes']

    static List<String> stable(Object table) {
        def lines = new File(table.toString()).readLines()
        def header = lines[0].split('\t', -1).toList()
        def keep = header.findIndexValues { column -> column in STABLE_COLUMNS }
        lines.collect { line ->
            def fields = line.split('\t', -1)
            keep.collect { index -> fields[index as int] }.join('\t')
        }
    }

    // The ΔPSI of each intron, by intron id, from the last column of an effect sizes table
    static Map<String, BigDecimal> deltapsi(Object table) {
        new File(table.toString()).readLines().drop(1).collectEntries { line ->
            def fields = line.split('\t', -1)
            [ fields[0], new BigDecimal(fields[-1]) ]
        }
    }

}
