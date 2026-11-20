import java.math.MathContext
import java.security.MessageDigest

/**
 * Reads CSV tables of statistics for snapshots that do not depend on the build.
 *
 * Builds of the same R packages differ in the last digits, so conda and the
 * containers disagree on the raw md5 of such tables. Decimal values are rounded
 * to a few significant digits, while integers and text are kept as they are.
 */
class RoundedCsv {

    static List<String> read(Object csv, int digits = 6) {
        def context = new MathContext(digits)
        new File(csv.toString()).readLines().collect { line ->
            line.split(',', -1).collect { value ->
                if (value ==~ /-?[0-9]+/ || !value.isBigDecimal()) {
                    return value
                }
                new BigDecimal(value).round(context).toString()
            }.join(',')
        }
    }

    static String md5(Object csv, int digits = 6) {
        def bytes = read(csv, digits).join('\n').getBytes('UTF-8')
        MessageDigest.getInstance('MD5').digest(bytes).encodeHex().toString()
    }

}
