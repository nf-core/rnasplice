//
// LeafCutter intron clustering and differential splicing subworkflow
//

include { REGTOOLS_JUNCTIONSEXTRACT as REGTOOLS_JUNCTIONSEXTRACT_FR } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { REGTOOLS_JUNCTIONSEXTRACT as REGTOOLS_JUNCTIONSEXTRACT_RF } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { REGTOOLS_JUNCTIONSEXTRACT as REGTOOLS_JUNCTIONSEXTRACT_XS } from '../../../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTERREGTOOLS } from '../../../modules/nf-core/leafcutter/clusterregtools/main'
include { LEAFCUTTER_DIFFERENTIALSPLICING } from '../../../modules/nf-core/leafcutter/differentialsplicing/main'
include { LEAFCUTTER_GTFTOEXONS } from '../../../modules/local/leafcutter/gtftoexons/main'
include { STRANDJUNCTIONS } from '../../../modules/local/strandjunctions/main'
include { resolveBlockingVariables } from '../utils_nfcore_rnasplice_pipeline'

workflow LEAFCUTTER {
    take:
    ch_genome_bam // channel: [ val(meta), path(bam) ], meta.condition is required
    ch_genome_bam_index // channel: [ val(meta), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gtf // channel: [ val(meta), path(gtf) ]
    ch_contrastsheet // channel: [ contrast:, treatment:, control:, blocking: (when set) ]
    ch_samplesheet // channel: file(samplesheet), read for the values of the blocking variables
    min_samples_per_condition // integer: fewest samples each condition of a contrast needs for leafcutter-ds to test it

    main:

    //
    // MODULE: REGTOOLS_JUNCTIONSEXTRACT
    //

    // `leafcutter-cluster` discards every junction that has no strand, so regtools has to
    // be told where to take the strand from. A stranded library has it in the orientation
    // of its reads, which regtools reads off the BAM: `FR` for a forward (second strand)
    // library, `RF` for a reverse (first strand) one. An unstranded library has no
    // orientation to read it from, so it is run with `XS`, which takes the strand from the
    // aligner tag of that name where there is one, and `STRANDJUNCTIONS` then fills in
    // the junctions left strandless from the splice motif and the annotation
    ch_bam_bai = ch_genome_bam
        .join(ch_genome_bam_index)
        .branch { meta, _bam, _bai ->
            forward: meta.strandedness == 'forward'
            reverse: meta.strandedness == 'reverse'
            unstranded: true
        }

    REGTOOLS_JUNCTIONSEXTRACT_FR(ch_bam_bai.forward, 'FR')
    REGTOOLS_JUNCTIONSEXTRACT_RF(ch_bam_bai.reverse, 'RF')
    REGTOOLS_JUNCTIONSEXTRACT_XS(ch_bam_bai.unstranded, 'XS')

    //
    // MODULE: STRANDJUNCTIONS
    //

    STRANDJUNCTIONS(REGTOOLS_JUNCTIONSEXTRACT_XS.out.junc, ch_fasta, ch_gtf)

    ch_junc = REGTOOLS_JUNCTIONSEXTRACT_FR.out.junc
        .mix(REGTOOLS_JUNCTIONSEXTRACT_RF.out.junc)
        .mix(STRANDJUNCTIONS.out.junc)

    //
    // MODULE: LEAFCUTTER_CLUSTERREGTOOLS
    //

    // The order the junction files are collected in does not reach the results:
    // `LEAFCUTTER_CLUSTERREGTOOLS` sorts the file list itself before handing it to
    // `leafcutter-cluster`, so no sorting is needed here
    ch_juncs = ch_junc
        .map { _meta, junc -> junc }
        .collect()
        .map { junc_files -> [[id: 'lc'], junc_files] }

    LEAFCUTTER_CLUSTERREGTOOLS(ch_juncs)

    //
    // MODULE: LEAFCUTTER_GTFTOEXONS
    //

    // The exon table `leafcutter-ds` labels each cluster with the genes it overlaps from
    LEAFCUTTER_GTFTOEXONS(ch_gtf)

    // A value channel, so that every contrast gets the exon table even if ch_gtf is a queue channel
    ch_exons = LEAFCUTTER_GTFTOEXONS.out.exons.first()

    //
    // MODULE: LEAFCUTTER_DIFFERENTIALSPLICING
    //

    // The samples of each condition, by the name `leafcutter-cluster` gives their column in
    // the counts table, which is the sample id the junction file is named after
    ch_condition_samples = ch_genome_bam
        .map { meta, _bam -> [meta.condition, meta.id] }
        .groupTuple(sort: true)
        .toList()
        .map { conditions -> [conditions.collectEntries()] }

    // The samplesheet row of each sample, the first for a sample with several runs
    ch_sample_rows = ch_samplesheet
        .splitCsv(header: true)
        .unique { row -> row.sample }
        .toList()
        .map { rows -> [rows.collectEntries { row -> [row.sample, row] }] }

    // `leafcutter-ds` stops with an error when either condition has fewer samples than
    // it needs to test a cluster. That is common in small designs, which can still be
    // clustered, so such a contrast is left out with a warning rather than failing the run
    ch_ds_contrasts = ch_contrastsheet
        .combine(ch_condition_samples)
        .combine(ch_sample_rows)
        .map { row, samples, rows_by_sample ->
            [row.treatment, row.control].each { condition ->
                if (!samples.containsKey(condition)) {
                    error("LeafCutter contrast '${row.contrast}': condition '${condition}' has no sample in the samplesheet")
                }
            }
            def groups = [row.treatment, row.control].collectEntries { condition -> [condition, samples[condition]] }
            // The blocking variables go into the groups file as extra columns, which
            // leafcutter-ds takes as confounders
            def blocking = []
            if (row.blocking) {
                def sample_rows = groups.collectMany { condition, ids -> ids.collect { id -> rows_by_sample[id] + [condition: condition] } }
                blocking = resolveBlockingVariables(row.contrast, 'LeafCutter', row.blocking, sample_rows)
            }
            def meta = [id: row.contrast, treatment: row.treatment, control: row.control]
            // Prefixed with the column name, as leafcutter-ds reads numeric values as continuous
            def confounders = groups.values().flatten().collectEntries { id ->
                [id, blocking.collect { column -> "${column}_${rows_by_sample[id][column].toString().trim()}".replaceAll(/\s+/, '_') }]
            }
            return [meta, groups, confounders]
        }
        .filter { meta, groups, _confounders ->
            def too_small = groups.findAll { _condition, ids -> ids.size() < min_samples_per_condition }
            if (too_small) {
                def counts = too_small.collect { condition, ids -> "'${condition}' has ${ids.size()}" }.join(', ')
                log.warn("LeafCutter contrast '${meta.id}' is not tested for differential splicing: each condition needs at least ${min_samples_per_condition} samples, and ${counts}. Lower --leafcutter_ds_min_samples_per_intron and --leafcutter_ds_min_samples_per_group to test it.")
            }
            return !too_small
        }

    // The groups file of each contrast, one sample per line followed by its condition and
    // its blocking variables
    ch_groups = ch_ds_contrasts
        .collectFile { meta, groups, confounders ->
            def lines = groups.collectMany { condition, ids -> ids.collect { id -> ([id, condition] + confounders[id]).join('\t') + '\n' } }
            ["${meta.id}_groups.txt", lines.join('')]
        }
        .map { groups_file -> [groups_file.name.replaceFirst(/_groups\.txt$/, ''), groups_file] }

    // `leafcutter-ds` reads the counts of each intron, which is the numerators table
    ch_ds_input = ch_ds_contrasts
        .map { meta, _groups, _confounders -> [meta.id, meta] }
        .join(ch_groups, failOnDuplicate: true, failOnMismatch: true)
        .combine(LEAFCUTTER_CLUSTERREGTOOLS.out.numers)
        .map { _id, meta, groups_file, _meta_lc, numers -> [meta, numers, groups_file] }

    LEAFCUTTER_DIFFERENTIALSPLICING(ch_ds_input, ch_exons)

    emit:
    juncs = ch_junc // channel: [ val(meta), path(junc) ]
    strand_summary = STRANDJUNCTIONS.out.summary // channel: [ val(meta), path(tsv) ], unstranded samples only
    counts = LEAFCUTTER_CLUSTERREGTOOLS.out.counts // channel: [ val(meta), path(counts) ]
    numers = LEAFCUTTER_CLUSTERREGTOOLS.out.numers // channel: [ val(meta), path(numers) ]
    groups = ch_groups // channel: [ val(contrast), path(txt) ], the groups file of each tested contrast, with its blocking variables
    pooled = LEAFCUTTER_CLUSTERREGTOOLS.out.pooled // channel: [ val(meta), path(pooled) ]
    refined = LEAFCUTTER_CLUSTERREGTOOLS.out.refined // channel: [ val(meta), path(refined) ]
    sortedlibs = LEAFCUTTER_CLUSTERREGTOOLS.out.sortedlibs // channel: [ val(meta), path(sortedlibs) ]
    exons = LEAFCUTTER_GTFTOEXONS.out.exons // channel: [ val(meta), path(exons) ]
    cluster_significance = LEAFCUTTER_DIFFERENTIALSPLICING.out.cluster_significance // channel: [ val(meta), path(txt) ]
    effect_sizes = LEAFCUTTER_DIFFERENTIALSPLICING.out.effect_sizes // channel: [ val(meta), path(txt) ]
}
