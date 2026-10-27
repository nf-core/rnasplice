//
// ASpli differential bin and junction usage subworkflow
//

include { ASPLI_COUNTS } from '../../../modules/local/aspli/counts'
include { ASPLI_DU     } from '../../../modules/local/aspli/du'

workflow ASPLI {
    take:
    ch_bam // channel: [ val(meta), path(bam) ], meta.condition is required
    ch_bai // channel: [ val(meta), path(bai) ]
    ch_gtf // channel: [ val(meta), path(gtf) ]
    ch_contrastsheet // channel: [ contrast:, treatment:, control: ]
    read_length // integer: --minReadLength of ASpli, 0 to estimate it from the BAM files

    main:

    //
    // The samples of the analysis, sorted by sample id so that the task hash does not
    // depend on the order the alignments finished in
    //
    ch_samples = ch_bam
        .join(ch_bai)
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .map { rows ->
            def conditions = rows.collect { meta, _bam, _bai -> meta.condition }.unique()
            if (conditions.size() < 2) {
                error("ASpli needs samples from at least two conditions, the samplesheet has: ${conditions.join(', ')}")
            }
            // ASpli counts every sample with one library type
            if (rows.collect { meta, _bam, _bai -> meta.single_end }.unique().size() > 1) {
                error("Please check input samplesheet -> Cannot run ASpli with mixed single and paired end samples.")
            }
            return rows
        }

    // The targets table ASpli reads the samples from, with the BAM files by the name
    // they are staged under
    ch_targets = ch_samples
        .map { rows ->
            def header = ['sample', 'condition', 'bam', 'single_end', 'strandedness'].join('\t')
            def lines = rows.collect { meta, bam, _bai -> [meta.id, meta.condition, bam.name, meta.single_end, meta.strandedness].join('\t') }
            return ([header] + lines).join('\n') + '\n'
        }
        .collectFile(name: 'aspli_targets.tsv')

    //
    // MODULE: ASPLI_COUNTS
    //

    // The genome is binned and the reads counted once, for all the contrasts
    ch_counts_input = ch_samples
        .combine(ch_targets)
        .map { items ->
            def targets = items.last()
            def rows = items.init()
            [[id: 'aspli'], targets, rows.collect { _meta, bam, _bai -> bam }, rows.collect { _meta, _bam, bai -> bai }]
        }

    ASPLI_COUNTS(
        ch_counts_input,
        ch_gtf,
        read_length,
    )

    //
    // MODULE: ASPLI_DU
    //

    // Every condition of a contrast needs samples, which is checked as soon as the
    // samples are known rather than once the reads are counted
    ch_conditions = ch_bam
        .map { meta, _bam -> meta.condition }
        .unique()
        .toSortedList()
        .map { conditions -> [conditions: conditions] }

    ch_contrasts = ch_contrastsheet
        .combine(ch_conditions)
        .map { row, samples ->
            [row.treatment, row.control].each { condition ->
                if (!samples.conditions.contains(condition)) {
                    error("ASpli contrast '${row.contrast}': condition '${condition}' has no sample in the samplesheet")
                }
            }
            return [id: row.contrast, treatment: row.treatment, control: row.control]
        }

    ch_du_input = ASPLI_COUNTS.out.features
        .join(ASPLI_COUNTS.out.gbcounts)
        .join(ASPLI_COUNTS.out.jcounts)
        .combine(ch_contrasts)
        .map { _meta, features, gbcounts, jcounts, contrast -> [contrast, features, gbcounts, jcounts] }

    ASPLI_DU(ch_du_input)

    emit:
    features = ASPLI_COUNTS.out.features // channel: [ val(meta), path(rds) ]
    gbcounts = ASPLI_COUNTS.out.gbcounts // channel: [ val(meta), path(rds) ]
    jcounts = ASPLI_COUNTS.out.jcounts // channel: [ val(meta), path(rds) ]
    counts = ASPLI_COUNTS.out.counts // channel: [ val(meta), [ path(tsv) ] ]
    junction_usage = ASPLI_COUNTS.out.junction_usage // channel: [ val(meta), [ path(tsv) ] ]
    settings = ASPLI_COUNTS.out.settings // channel: [ val(meta), path(tsv) ]
    gb_du = ASPLI_DU.out.gb_du // channel: [ val(meta), [ path(tsv) ] ]
    j_du = ASPLI_DU.out.j_du // channel: [ val(meta), [ path(tsv) ] ]
    splicing_report = ASPLI_DU.out.splicing_report // channel: [ val(meta), [ path(tsv) ] ]
    integrated_signals = ASPLI_DU.out.integrated_signals // channel: [ val(meta), path(tsv) ]
    integrated_signals_html = ASPLI_DU.out.integrated_signals_html // channel: [ val(meta), path(dir) ]
    rds = ASPLI_DU.out.rds // channel: [ val(meta), path(rds) ]
}
