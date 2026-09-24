//
// rMATS differential splicing analysis
//

include { RMATS_PREP } from '../../../modules/nf-core/rmats/prep'

include { CREATE_BAMLIST } from '../../../modules/local/create_bamlist'
include { RMATS_POST } from '../../../modules/local/rmats_post'

workflow RMATS {
    take:
    ch_samples // channel: [ val(sample_id), val(condition) ], in samplesheet order
    ch_contrastsheet // channel: [ contrast:, treatment:, control: ]
    ch_bam // channel: [ val(meta), path(bam) ], meta.condition is required
    ch_gtf // channel: [ val(meta), path(gtf) ]
    rmats_read_len // integer: --readLength of both rMATS steps
    rmats_paired_stats // boolean: pair the treatment and control samples

    main:

    //
    // MODULE: RMATS_PREP
    //

    // The prep step reads each BAM file once, whatever the number of contrasts it takes
    // part in, and records it in a `.rmats` file under the BAM file name
    RMATS_PREP(
        ch_bam,
        ch_gtf,
        rmats_read_len,
    )

    //
    // Sample ids grouped by condition, in samplesheet order, which is what the paired
    // model pairs treatment and control samples by. This is the only view of the
    // samplesheet rMATS needs: how many conditions there are, and the order of the
    // samples within each of them.
    //
    ch_condition_samples = ch_samples
        .toList()
        .map { rows ->
            def by_condition = [:]
            rows.each { sample_id, condition ->
                def ids = by_condition.computeIfAbsent(condition) { [] }
                if (!ids.contains(sample_id)) {
                    ids << sample_id
                }
            }
            return by_condition
        }

    //
    // The samples of each rMATS run: [ contrast, cond1, ids1, cond2, ids2 ]
    //

    // With one condition there is nothing to contrast, so it is profiled on its own and
    // the second condition stays empty
    ch_single_samples = ch_condition_samples.flatMap { by_condition ->
        by_condition.size() == 1
            ? by_condition.collect { condition, ids -> ["${condition}_profiling".toString(), condition, ids, '', []] }
            : []
    }

    // The contrastsheet only makes sense with more than one condition, so its rows are
    // dropped in single condition mode
    ch_contrast_samples = ch_contrastsheet
        .combine(ch_condition_samples)
        .filter { _row, by_condition -> by_condition.size() > 1 }
        .map { row, by_condition ->
            [row.treatment, row.control].each { condition ->
                if (!by_condition.containsKey(condition)) {
                    error("rMATS contrast '${row.contrast}': condition '${condition}' has no sample in the samplesheet")
                }
            }
            def ids1 = by_condition[row.treatment]
            def ids2 = by_condition[row.control]
            // The paired model matches the n-th treatment sample with the n-th control one
            if (rmats_paired_stats && ids1.size() != ids2.size()) {
                error("Paired rMATS contrast '${row.contrast}': unequal sample counts, ${row.treatment}: ${ids1}, ${row.control}: ${ids2}. Each treatment sample needs a control counterpart.")
            }
            return [row.contrast, row.treatment, ids1, row.control, ids2]
        }

    // Exactly one of the two is populated, so rMATS runs off a single set of tuples
    // whichever mode the samplesheet puts it in
    ch_run_samples = ch_single_samples.mix(ch_contrast_samples)

    // One entry per sample of a run: [ sample_id, contrast, condition number, position ]
    ch_run_sample_index = ch_run_samples.flatMap { contrast, _cond1, ids1, _cond2, ids2 ->
        def entries1 = ids1.withIndex().collect { sample_id, idx -> [sample_id, contrast, 1, idx] }
        def entries2 = ids2.withIndex().collect { sample_id, idx -> [sample_id, contrast, 2, idx] }
        return entries1 + entries2
    }

    //
    // MODULE: CREATE_BAMLIST
    //

    // The BAM files of each run, in samplesheet order within each condition
    ch_run_bams = ch_run_sample_index
        .combine(ch_bam.map { meta, bam -> [meta.id, bam] }, by: 0)
        .map { _sample_id, contrast, cond_num, idx, bam -> [contrast, [cond_num, idx, bam]] }
        .groupTuple()
        .map { contrast, entries ->
            def sorted = entries.sort { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }
            def bam1 = sorted.findAll { entry -> entry[0] == 1 }.collect { entry -> entry[2] }
            def bam2 = sorted.findAll { entry -> entry[0] == 2 }.collect { entry -> entry[2] }
            [contrast, bam1, bam2]
        }

    ch_bamlist_input = ch_run_samples
        .join(ch_run_bams, by: 0)
        .map { contrast, cond1, _ids1, cond2, _ids2, bam1, bam2 ->
            [contrast, cond1, bam1, cond2, bam2]
        }

    CREATE_BAMLIST(
        ch_bamlist_input
    )

    //
    // MODULE: RMATS_POST
    //

    // The post step of a run takes the `.rmats` files of every sample in its bam lists.
    // The bam lists carry the BAM file names, which is what rMATS matches the `.rmats`
    // files by, so the samples of a run are looked up by sample id here
    ch_run_rmats = ch_run_sample_index
        .map { sample_id, contrast, _cond_num, _idx -> [sample_id, contrast] }
        .combine(RMATS_PREP.out.rmats.map { meta, rmats -> [meta.id, rmats] }, by: 0)
        .map { _sample_id, contrast, rmats -> [contrast, rmats] }
        .groupTuple()

    // CREATE_BAMLIST only writes the second bam list when there is a second condition,
    // so a single condition run joins with no bam list 2. It has to become an empty list
    // rather than an empty string, which a `path` input rejects
    ch_post_input = ch_run_samples
        .join(CREATE_BAMLIST.out.bam_list1, by: 0)
        .join(CREATE_BAMLIST.out.bam_list2, by: 0, remainder: true)
        .join(ch_run_rmats, by: 0)
        .map { contrast, cond1, _ids1, cond2, _ids2, bam_list1, bam_list2, rmats ->
            // A single condition run has no control, so its meta map has no control key
            def meta = cond2 ? [id: contrast, treatment: cond1, control: cond2] : [id: contrast, treatment: cond1]
            // Sorted so that the task hash does not depend on the order the prep tasks
            // finished in
            [meta, rmats.sort { rmats_file -> rmats_file.name }, bam_list1, bam_list2 ?: []]
        }

    RMATS_POST(
        ch_post_input,
        ch_gtf,
        rmats_read_len,
    )

    emit:
    rmats = RMATS_PREP.out.rmats // channel: [ val(meta), path(rmats) ]
    read_outcomes = RMATS_PREP.out.read_outcomes // channel: [ val(meta), path(txt) ]
    bam_list1 = CREATE_BAMLIST.out.bam_list1 // channel: [ val(contrast), path(txt) ]
    bam_list2 = CREATE_BAMLIST.out.bam_list2 // channel: [ val(contrast), path(txt) ]
    mats = RMATS_POST.out.mats // channel: [ val(meta), [ path(txt) ] ]
    from_gtf = RMATS_POST.out.from_gtf // channel: [ val(meta), [ path(txt) ] ]
    raw_input = RMATS_POST.out.raw_input // channel: [ val(meta), [ path(txt) ] ]
    summary = RMATS_POST.out.summary // channel: [ val(meta), path(txt) ]
    post_log = RMATS_POST.out.log // channel: [ val(meta), path(log) ]
}
