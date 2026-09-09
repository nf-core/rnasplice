//
// rMATS differential splicing analysis
//

include { CREATE_BAMLIST         } from '../../../modules/local/create_bamlist'
include { RMATS_PREP             } from '../../../modules/local/rmats_prep'
include { RMATS_POST             } from '../../../modules/local/rmats_post'


workflow RMATS {
    take:
    ch_samples               // channel: [ sample_id, condition ], in samplesheet order
    ch_contrastsheet         // channel: [ contrast, treatment, control ]
    ch_genome_bam_conditions // channel: [ condition, meta, bam ]
    gtf                      // channel: path(genome.gtf)
    rmats_read_len
    rmats_splice_diff_cutoff
    rmats_novel_splice_site
    rmats_min_intron_len
    rmats_max_exon_len
    rmats_paired_stats

    main:

    //
    // Samples grouped by condition, keeping the samplesheet order, which is what the
    // paired model pairs treatment and control samples by. Technical replicates share
    // a sample id, so they collapse into a single entry.
    //
    // This is the only view of the samplesheet rMATS needs: how many conditions there
    // are, and the order of the samples within each of them.
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
    // SINGLE CONDITION MODE
    //
    // With one condition there is nothing to contrast, so it is profiled on its own and
    // the second half of the tuple stays empty.
    //
    ch_single_bamlist = ch_genome_bam_conditions
        .groupTuple(by: 0)
        .combine( ch_condition_samples )
        .filter { _condition, _metas, _bams, by_condition -> by_condition.size() == 1 }
        .map { condition, metas, bams, _by_condition ->
            [ "${condition}_profiling", condition, metas, bams, '', [], [] ]
        }

    //
    // TWO CONDITIONS MODE
    //
    // The contrastsheet only makes sense with more than one condition, so its rows are
    // dropped in single condition mode.
    //
    ch_contrasts = ch_contrastsheet
        .combine( ch_condition_samples )
        .filter { _row, by_condition -> by_condition.size() > 1 }
        .map { row, _by_condition -> row }

    if (rmats_paired_stats) {

        //
        // PAIRED SAMPLES
        //

        // Position of each sample within its condition, which is how a treatment sample
        // is matched to its control counterpart
        ch_sample_index = ch_condition_samples
            .flatMap { by_condition ->
                by_condition.collectMany { _condition, sample_ids ->
                    sample_ids.withIndex().collect { sample_id, idx -> [ sample_id, idx ] }
                }
            }

        ch_genome_bam_conditions
            .multiMap { condition, meta, bam ->
                tx:   [ condition, meta, bam ]
                ctrl: [ condition, meta, bam ]
            }
            .set { ch_bams_fork }

        ch_indexed_bams_tx = ch_bams_fork.tx
            .map { condition, meta, bam -> [ meta.id, condition, meta, bam ] }
            .join( ch_sample_index, by: 0 )
            .map { _id, condition, meta, bam, idx -> [ condition, idx, meta, bam ] }

        ch_indexed_bams_ctrl = ch_bams_fork.ctrl
            .map { condition, meta, bam -> [ meta.id, condition, meta, bam ] }
            .join( ch_sample_index, by: 0 )
            .map { _id, condition, meta, bam, idx -> [ condition, idx, meta, bam ] }


        // Build per-contrast pairs keyed by contrast__sample_id
        ch_tx = ch_contrasts
            .map { row -> [ row.treatment, row ] }           // key: treatment condition
            .combine( ch_indexed_bams_tx, by: 0 )            // match BAMs by condition
            .map { _condition, row, idx, meta, bam ->
                [ row.contrast, meta.id, idx, row, meta, bam ]  // key: contrast + sample id (separate)
            }

        ch_ctrl = ch_contrasts
            .map { row -> [ row.control, row ] }             // key: control condition
            .combine( ch_indexed_bams_ctrl, by: 0 )          // match BAMs by condition
            .map { _condition, row, idx, meta, bam ->
                [ row.contrast, meta.id, idx, row, meta, bam ]  // same contrast + sample id
            }

        // Join treatment and control on contrast__sample_id
        ch_fully_paired = ch_tx
            .join( ch_ctrl, by: [0, 2] )  // join on both contrast AND sample id
            .map { _contrast, _sample_id, tx_idx, tx_row, tx_meta, tx_bam, ctrl_idx, _ctrl_row, ctrl_meta, ctrl_bam ->
                tx_row + [
                    tx_idx   : tx_idx,
                    tx_meta  : tx_meta,
                    tx_bam   : tx_bam,
                    ctrl_idx : ctrl_idx,
                    ctrl_meta: ctrl_meta,
                    ctrl_bam : ctrl_bam
                ]
            }

        ch_contrasts_bamlist = ch_fully_paired
            .map { it -> [ it.contrast, it ] }
            .groupTuple(by: 0)
            .map { contrast, pairs ->
                def cond1  = pairs[0].treatment
                def cond2  = pairs[0].control
                def sorted = pairs.sort { it -> it.tx_idx }
                def meta1  = sorted.collect { it -> it.tx_meta }
                def bam1   = sorted.collect { it -> it.tx_bam }
                def meta2  = sorted.collect { it -> it.ctrl_meta }
                def bam2   = sorted.collect { it -> it.ctrl_bam }

                if (meta1.size() != meta2.size()) {
                    error(
                        "Paired rMATS contrast '${contrast}': unequal sample counts — " +
                        "${cond1}: ${meta1.collect { it -> it.id }}, " +
                        "${cond2}: ${meta2.collect { it -> it.id }}. " +
                        "Each sample ID must appear in both conditions."
                    )
                }

                return [ contrast, cond1, meta1, bam1, cond2, meta2, bam2 ]
            }

    } else {

        //
        // UNPAIRED SAMPLES
        //

        ch_grouped_bams = ch_genome_bam_conditions
            .groupTuple(by: 0)

        ch_contrasts_bamlist = ch_contrasts
            .map { row -> [ row.treatment, row ] }
            .join( ch_grouped_bams, by: 0 )
            .map { _tx, row, tx_metas, tx_bams ->
                [ row.control, row + [ tx_metas: tx_metas, tx_bams: tx_bams ] ]
            }
            .join( ch_grouped_bams, by: 0 )
            .map { _ctrl, row, ctrl_metas, ctrl_bams ->
                [ row.contrast, row.treatment, row.tx_metas, row.tx_bams, row.control, ctrl_metas, ctrl_bams ]
            }
    }

    // Exactly one of the two is populated, so rMATS runs off a single set of tuples
    // whichever mode the samplesheet puts it in
    ch_all_contrasts_bamlist = ch_single_bamlist.mix( ch_contrasts_bamlist )

    CREATE_BAMLIST(
        ch_all_contrasts_bamlist
            .map { contrast, cond1, _meta1, bam1, cond2, _meta2, bam2 ->
                [ contrast, cond1, bam1, cond2, bam2 ]
            }
    )

    // CREATE_BAMLIST only writes the second bam list when there is a second condition,
    // so single condition contrasts join with no bam list 2. It has to become an empty
    // list rather than an empty string, which a `path` input rejects
    ch_prep_ready = ch_all_contrasts_bamlist
        .join( CREATE_BAMLIST.out.bam_list1, by: 0 )
        .join( CREATE_BAMLIST.out.bam_list2, by: 0, remainder: true )
        .map { contrast, cond1, meta1, bam1, cond2, meta2, bam2, bam1_txt, bam2_txt ->
            return [ contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt ?: [] ]
        }

    RMATS_PREP(
        gtf,
        ch_prep_ready,
        rmats_read_len,
        rmats_splice_diff_cutoff,
        rmats_novel_splice_site,
        rmats_min_intron_len,
        rmats_max_exon_len,
    )

    ch_post_ready = ch_prep_ready
        .join( RMATS_PREP.out.rmats_temp, by: 0 )
        .map { contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt, rmats_temp ->
            [ contrast, cond1, meta1, bam1, bam1_txt, cond2, meta2, bam2, bam2_txt, rmats_temp ]
        }

    RMATS_POST(
        gtf,
        ch_post_ready,
        rmats_read_len,
        rmats_splice_diff_cutoff,
        rmats_novel_splice_site,
        rmats_min_intron_len,
        rmats_max_exon_len,
        rmats_paired_stats,
    )

    ch_rmats_prep     = RMATS_PREP.out.rmats_temp
    ch_rmats_prep_log = RMATS_PREP.out.log
    ch_rmats_post     = RMATS_POST.out.rmats_post
    ch_rmats_post_log = RMATS_POST.out.log

    emit:
    rmats_prep     = ch_rmats_prep
    rmats_prep_log = ch_rmats_prep_log
    rmats_post     = ch_rmats_post
    rmats_post_log = ch_rmats_post_log
}
