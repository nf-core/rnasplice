//
// Subworkflow with functionality specific to the nf-core/rnasplice pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    source            //  string: Type of input data [fastq, genome_bam, transcriptome_bam, salmon_results]
    contrasts         //  string: Path to the contrasts sheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/rnasplice ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/rnasplice/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    // The samplesheet columns depend on --source, so it is validated against the schema
    // matching the declared source, selected by `samplesheetSchema`, which errors on an
    // unknown --source. This is the only place the pipeline itself reads the samplesheet:
    // what is emitted here is a parsed channel for the consumers that need the values and
    // the validated path for the processes that read the sheet themselves.
    //

    def samplesheet_rows = samplesheetToList(input, samplesheetSchema(source))

    if (source == 'fastq') {
        channel
            .fromList(samplesheet_rows)
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }
    } else if (source == 'transcriptome_bam') {
        // The samplesheet carries both genome_bam and transcriptome_bam; only the
        // transcriptome BAM is read from here. Strandedness and read type cannot be
        // derived from a BAM file, so they come from the samplesheet, which defaults
        // them to 'unstranded' and paired end when the columns are missing
        channel
            .fromList(samplesheet_rows)
            .map {
                meta, _genome_bam, transcriptome_bam ->
                    return [ [ id: meta.id, condition: meta.condition, strandedness: meta.strandedness, single_end: meta.single_end ], [ transcriptome_bam ] ]
            }
            .set { ch_samplesheet }
    } else if (source == 'genome_bam') {
        channel
            .fromList(samplesheet_rows)
            .map {
                meta, genome_bam ->
                    return [ [ id: meta.id, condition: meta.condition, strandedness: meta.strandedness, single_end: meta.single_end ], [ genome_bam ] ]
            }
            .set { ch_samplesheet }
    } else {
        // salmon_results carries a single path column, and the quantification it holds
        // has already been done, so no strandedness is needed
        channel
            .fromList(samplesheet_rows)
            .map {
                meta, input_file ->
                    return [ [ id: meta.id, condition: meta.condition ], [ input_file ] ]
            }
            .set { ch_samplesheet }
    }

    //
    // Create channel from the contrasts file provided through params.contrasts
    //
    // Validated here next to the samplesheet, so each sheet is read once and both are
    // emitted the same way: the rows for the consumers that work from the values, the
    // validated path for the processes that read the sheet themselves.
    //

    channel
        .fromList(samplesheetToList(contrasts, "${projectDir}/assets/schema_contrasts.json"))
        .map { row ->
            // Every column of the contrastsheet is a meta field, so a row comes back as a
            // single element list holding the meta map, not as the map itself
            def meta = row[0]
            validateInputContrastsheet([[meta]])
            return [ contrast: meta.contrast, treatment: meta.treatment, control: meta.control ]
        }
        .set { ch_contrastsheet }

    emit:
    samplesheet        = ch_samplesheet                                      // channel: [ val(meta), [ files ] ]
    samplesheet_file   = channel.value(file(input, checkIfExists: true))     // channel: path(samplesheet.csv)
    contrastsheet      = ch_contrastsheet                                    // channel: [ contrast, treatment, control ]
    contrastsheet_file = channel.value(file(contrasts, checkIfExists: true)) // channel: path(contrastsheet.csv)
    versions           = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    if (!params.fasta) {
        error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
    }

    if (!params.gtf && !params.gff) {
        error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
    }

    if (params.gtf) {
        if (params.gff) {
            gtfGffWarn()
        }
        if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
            ncbiGenomeWarn()
        }
        if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
            ucscGenomeWarn()
        }
    }

    if (params.transcript_fasta) {
        transcriptsFastaWarn()
    }

    def valid_params = [aligners: ['star', 'star_salmon'], pseudoaligners: ['salmon']]

    if (!params.skip_alignment) {
        if (!valid_params['aligners'].contains(params.aligner)) {
            error("Invalid option: '${params.aligner}'. Valid options for '--aligner': '${valid_params['aligners'].join(', ')}'.")
        }
    }
    else {
        if (!params.pseudo_aligner) {
            error("--skip_alignment specified without --pseudo_aligner...please specify e.g. --pseudo_aligner '${valid_params['pseudoaligners'][0]}'.")
        }
        skipAlignmentWarn()
    }

    if (params.pseudo_aligner) {
        if (!valid_params['pseudoaligners'].contains(params.pseudo_aligner)) {
            error("Invalid option: '${params.pseudo_aligner}'. Valid options for '--pseudo_aligner': '${valid_params['pseudoaligners'].join(', ')}'.")
        }
        else {
            if (!(params.salmon_index || params.transcript_fasta || (params.fasta && (params.gtf || params.gff)))) {
                error("To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff.")
            }
        }
    }

    if (params.clusterevents_local_event && !params.diffsplice_local_event) {
        error("--clusterevents_local_event specified without --diffsplice_local_event... please specify e.g. --diffsplice_local_event=true")
    }

    if (params.clusterevents_isoform && !params.diffsplice_isoform) {
        error("--clusterevents_isoform specified without diffsplice_isoform... please specify e.g. --diffsplice_isoform=true")
    }

    if (!params.skip_alignment) {
        if (params.aligner == "star_salmon" && params.pseudo_aligner == "salmon") {
            log.warn("Both --aligner=star_salmon and --pseudo_aligner=salmon specified. Downstream analyses will be performed on both salmon output files.")
        }
    }

    if (!params.input) {
        error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
    }
}

//
// Schema the input samplesheet is validated against, selected from --source
//
// Each source declares a different set of columns, so there is one schema per source.
// Every place that reads the samplesheet goes through this function, so a samplesheet
// can never be accepted by one consumer and rejected by another.
//
def samplesheetSchema(source) {
    def source_schemas = [
        fastq:             "${projectDir}/assets/schema_input.json",
        genome_bam:        "${projectDir}/assets/schema_input_genome_bam.json",
        transcriptome_bam: "${projectDir}/assets/schema_input_transcriptome_bam.json",
        salmon_results:    "${projectDir}/assets/schema_input_salmon_results.json",
    ]

    if (!source_schemas.containsKey(source)) {
        error("Invalid --source parameter: '${source}'. Must be one of: ${source_schemas.keySet().join(', ')}")
    }

    return source_schemas[source]
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastqs]
}

//
// Validate channels from input contrastsheet
//
def validateInputContrastsheet(input) {
    def metas = input[0]
    def meta = metas[0]

    if (meta.treatment == meta.control) {
        error("Please check input contrastsheet -> Treatment and control conditions cannot be the same: ${meta.contrast} (${meta.treatment} vs ${meta.control})")
    }

    return meta
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = ["Tools used in the workflow included:", "FastQC (Andrews 2010),", "MultiQC (Ewels et al. 2016)", "."].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = ["<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>", "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    def methods_text = mqc_methods_yaml.text
    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Function to check whether biotype field exists in GTF file
//
def biotypeInGtf(gtf_file, biotype) {
    def hits = 0
    gtf_file.eachLine { line ->
        def attributes = line.split('\t')[-1].split()
        if (attributes.contains(biotype)) {
            hits += 1
        }
    }
    if (hits) {
        return true
    }
    else {
        log.warn(
            "=============================================================================\n" + "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" + "  Biotype QC will be skipped to circumvent the issue below:\n" + "  https://github.com/nf-core/rnaseq/issues/460\n\n" + "  Amend '--featurecounts_group_type' to change this behaviour.\n" + "==================================================================================="
        )
        return false
    }
}

//
// Function to generate an error if contigs in genome fasta file > 512 Mbp
//
def checkMaxContigSize(fai_file) {
    def max_size = 512000000
    fai_file.eachLine { line ->
        def lspl = line.split('\t')
        def chrom = lspl[0]
        def size = lspl[1]
        if (size.toInteger() > max_size) {
            error(
                "=============================================================================\n" + "  Contig longer than ${max_size}bp found in reference genome!\n\n" + "  ${chrom}: ${size}\n\n" + "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" + "  Please see:\n" + "  https://github.com/nf-core/rnaseq/issues/744\n" + "============================================================================="
            )
        }
    }
}

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}

//
// Print a warning if using GRCh38 assembly from igenomes.config
//
def ncbiGenomeWarn() {
    log.warn(
        "=============================================================================\n" + "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" + "  Biotype QC will be skipped to circumvent the issue below:\n" + "  https://github.com/nf-core/rnaseq/issues/460\n\n" + "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" + "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" + "==================================================================================="
    )
}

//
// Print a warning if using a UCSC assembly from igenomes.config
//
def ucscGenomeWarn() {
    log.warn(
        "=============================================================================\n" + "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" + "  Biotype QC will be skipped to circumvent the issue below:\n" + "  https://github.com/nf-core/rnaseq/issues/460\n\n" + "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" + "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" + "==================================================================================="
    )
}

//
// Print a warning if both GTF and GFF have been provided
//
def gtfGffWarn() {
    log.warn(
        "=============================================================================\n" + "  Both '--gtf' and '--gff' parameters have been provided.\n" + "  Using GTF file as priority.\n" + "==================================================================================="
    )
}

//
// Print a warning if using '--transcript_fasta'
//
def transcriptsFastaWarn() {
    log.warn(
        "=============================================================================\n" + "  '--transcript_fasta' parameter has been provided.\n" + "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" + "  Please see:\n" + "  https://github.com/nf-core/rnaseq/issues/753\n" + "==================================================================================="
    )
}

//
// Print a warning if --skip_alignment has been provided
//
def skipAlignmentWarn() {
    log.warn(
        "=============================================================================\n" + "  '--skip_alignment' parameter has been provided.\n" + "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" + "==================================================================================="
    )
}

//
// Print a warning if LeafCutter is run on unstranded BAM input, whose junctions carry a
// strand only if the BAM files happen to have XS tags
//
def leafcutterUnstrandedBamWarn(samples) {
    samples
        .filter { meta, _files -> meta.strandedness == 'unstranded' }
        .map { meta, _files -> meta.id }
        .collect()
        .map { ids ->
            if (!ids) {
                return
            }
            log.warn(
                "=============================================================================\n" + "  LeafCutter takes the junction strand of an unstranded sample from the aligner\n" + "  XS tag, which the pipeline only asks STAR for when it does the alignment itself.\n" + "  The junctions of a BAM file without the tag have no strand and are discarded,\n" + "  which leaves the clusters empty. Samples: ${ids.join(', ')}.\n\n" + "  Set 'strandedness' in the samplesheet if the library is stranded.\n" + "==================================================================================="
            )
        }
}

//
// Exit pipeline if rMATS requested with mixed single and paired end samples
//
def rmatsReadError(samples) {
    samples
        .map { meta, _files -> meta.single_end }
        .unique()
        .collect()
        .map {
            if (it.size() > 1) {
                error("Please check input samplesheet -> Cannot run rMats with mixed single and paired end samples.")
            }
        }
}

//
// Exit pipeline if rMATS requested with mixed stranded samples
//
def rmatsStrandednessError(samples) {
    samples
        .map { meta, _files -> meta.strandedness }
        .unique()
        .collect()
        .map {
            if (it.size() > 1) {
                error("Please check input samplesheet -> Cannot run rMats with mixed stranded samples.")
            }
        }
}
