//
// IsoformSwitchAnalyzeR subworkflow: one run per contrast, with the blocking variables of
// the contrast as covariates of the switch test
//

include { ISOFORMSWITCHANALYZER     } from '../../../modules/local/isoformswitchanalyzer'
include { resolveBlockingVariables  } from '../utils_nfcore_rnasplice_pipeline'

workflow ISOFORMSWITCHANALYZER_CONTRASTS {
    take:
    ch_salmon_results  // channel: [ val(meta), path(quant_dir) ], meta.id is the sample
    ch_samplesheet     // channel: file(samplesheet), read for the conditions and blocking variables
    ch_contrastsheet   // channel: [ contrast:, treatment:, control:, blocking: (when set) ]
    gtf                // channel: path(gtf)
    transcript_fasta   // channel: path(fasta)
    alpha              // value: FDR cutoff
    dIF                // value: minimum absolute difference in isoform fraction

    main:

    // The samplesheet row of each sample, the first for a sample with several runs
    ch_sample_rows = ch_samplesheet
        .splitCsv(header: true)
        .unique { row -> row.sample }
        .toList()
        .map { rows -> [rows.collectEntries { row -> [row.sample, row] }] }

    // One design table and contrast table per contrast, from the samples of its two
    // conditions. Blocking variables are prefixed with their column name, so that
    // numeric values are not read as continuous
    ch_contrast_inputs = ch_contrastsheet
        .combine(ch_salmon_results.toList())
        .combine(ch_sample_rows)
        .map { row, salmon, rows_by_sample ->
            def conditions = [row.treatment, row.control]
            def selected = salmon.findAll { meta, _dir -> rows_by_sample[meta.id]?.condition in conditions }
            def sample_rows = selected.collect { meta, _dir -> rows_by_sample[meta.id] }
            def blocking = row.blocking ? resolveBlockingVariables(row.contrast, 'IsoformSwitchAnalyzeR', row.blocking, sample_rows) : []
            def design = [(['sample', 'condition'] + blocking).join(',')] + sample_rows.collect { sample_row ->
                ([sample_row.sample, sample_row.condition] + blocking.collect { column -> "${column}_${sample_row[column].toString().trim()}".replaceAll(/[\s,]+/, '_') }).join(',')
            }
            def contrast = ['contrast,treatment,control', "${row.contrast},${row.treatment},${row.control}"]
            [row.contrast, selected.collect { _meta, dir -> dir }, design.join('\n') + '\n', contrast.join('\n') + '\n']
        }

    ch_designs = ch_contrast_inputs
        .collectFile { id, _dirs, design, _contrast -> ["${id}_design.csv", design] }
        .map { design_file -> [design_file.name.replaceFirst(/_design\.csv$/, ''), design_file] }

    ch_contrasts = ch_contrast_inputs
        .collectFile { id, _dirs, _design, contrast -> ["${id}_contrast.csv", contrast] }
        .map { contrast_file -> [contrast_file.name.replaceFirst(/_contrast\.csv$/, ''), contrast_file] }

    ch_isa_input = ch_contrast_inputs
        .map { id, dirs, _design, _contrast -> [id, dirs] }
        .join(ch_designs, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_contrasts, failOnDuplicate: true, failOnMismatch: true)
        .map { id, dirs, design_file, contrast_file -> [[id: id], dirs, design_file, contrast_file] }

    ISOFORMSWITCHANALYZER(ch_isa_input, gtf, transcript_fasta, alpha, dIF)

    emit:
    summary         = ISOFORMSWITCHANALYZER.out.isoformswitchanalyzer_summary          // channel: [ val(meta), path(csv) ]
    isoform_features = ISOFORMSWITCHANALYZER.out.isoformswitchanalyzer_isoformFeatures // channel: [ val(meta), path(csv) ]
    switchlist_rds  = ISOFORMSWITCHANALYZER.out.switchlist_rds                         // channel: [ val(meta), path(rds) ]
    results         = ISOFORMSWITCHANALYZER.out.results                                // channel: [ val(meta), path(dir) ]
}
