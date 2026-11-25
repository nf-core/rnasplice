//
// MultiQC report assembly for nf-core/rnasplice
//

include { CUSTOM_SPLICINGSUMMARY } from '../../../modules/local/custom/splicingsummary'
include { MULTIQC                } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from 'plugin/nf-core-utils'
include { methodsDescriptionText } from '../utils_nfcore_rnasplice_pipeline'

workflow MULTIQC_RNASPLICE {

    take:
    ch_multiqc_files             // channel: [ val(meta), path(file) or [ path(files) ] ], the outputs MultiQC modules read as they are
    ch_fail_trimming_tsv         // channel: val(tsv), the samples that failed the trimmed reads threshold
    ch_strand_check_tsv          // channel: val(tsv), the strand check of each sample
    ch_leafcutter_strand_summary // channel: [ val(meta), path(tsv) ], unstranded samples only
    ch_rmats_read_outcomes       // channel: [ val(meta), path(txt) ]
    ch_rmats_summary             // channel: [ val(meta), path(txt) ], meta.id is the contrast
    ch_splicing_results          // channel: [ val(label), path(result) ], label: [ kind, route, id ], see CUSTOM_SPLICINGSUMMARY
    ch_collated_versions         // channel: path(versions.yml)
    multiqc_config               // string: path to a custom MultiQC config, or null
    multiqc_logo                 // string: path to a custom logo, or null
    multiqc_methods_description  // string: path to a custom methods description, or null

    main:

    def ch_workflow_summary = channel
        .value(paramsSummaryMultiqc(paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")))
        .collectFile(name: 'workflow_summary_mqc.yaml')

    def ch_methods_description_template = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel
        .value(methodsDescriptionText(ch_methods_description_template))
        .collectFile(name: 'methods_description_mqc.yaml', sort: true)

    // One table row per unstranded sample, sorted by sample so the report does not
    // depend on the order the tasks finished in
    def ch_leafcutter_strands = ch_leafcutter_strand_summary
        .map { _meta, summary -> summary }
        .collectFile(name: 'leafcutter_junction_strands_mqc.tsv', keepHeader: true, skip: 1, sort: { summary -> summary.name })

    // Every contrast writes a `summary.txt`, and MultiQC names it after its file
    def ch_rmats = ch_rmats_read_outcomes
        .map { _meta, txt -> txt }
        .mix(ch_rmats_summary.collectFile { meta, summary -> ["${meta.id}_summary.txt", summary] })

    // Sorted so that a resumed run finds the cached task
    CUSTOM_SPLICINGSUMMARY(
        ch_splicing_results
            .toSortedList { a, b -> a[0].toString() <=> b[0].toString() ?: a[1].name <=> b[1].name }
            .filter { results -> results }
            .map { results -> [[id: 'splicing_summary'], results.collect { label, _result -> label }, results.collect { _label, result -> result }] }
    )

    def ch_files = ch_multiqc_files
        .map { _meta, files -> files }
        .mix(ch_fail_trimming_tsv.collectFile(name: 'fail_trimmed_samples_mqc.tsv'))
        .mix(ch_strand_check_tsv.collectFile(name: 'fail_strand_check_mqc.tsv'))
        .mix(ch_leafcutter_strands)
        .mix(ch_rmats)
        .mix(CUSTOM_SPLICINGSUMMARY.out.sample_metrics.map { _meta, tsv -> tsv })
        .mix(CUSTOM_SPLICINGSUMMARY.out.differential.map { _meta, tsv -> tsv })
        .mix(ch_workflow_summary)
        .mix(ch_methods_description)
        .mix(ch_collated_versions)

    MULTIQC(
        ch_files.flatten().collect().map { files ->
            [
                [id: 'rnasplice'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: [ path(multiqc_report.html) ]
    data   = MULTIQC.out.data.map { _meta, data -> data }                  // channel: path(multiqc_data)
}
