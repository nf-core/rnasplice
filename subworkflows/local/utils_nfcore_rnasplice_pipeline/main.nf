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
        null
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
    // Every input type is an optional column, and which ones are filled decides how each
    // method gets its input, see `getInputRoutes`. This is where the samplesheet is
    // validated: what is emitted is a parsed channel for the consumers that need the values
    // and a sample sheet for the processes that read one themselves.
    //

    def samplesheet_rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    validateInputRoutes(getInputRoutes(getSamplesheetInputTypes(samplesheet_rows)))

    // A sample may have one row per sequencing run. Its runs are collected under its id, and
    // the FASTQ files of all of them kept for concatenation, while the BAM files and Salmon
    // results, which describe the whole sample, must be the same on every row
    channel
        .fromList(samplesheet_rows)
        .map { meta, fastq_1, fastq_2, genome_bam, transcriptome_bam, salmon_results ->
            def fastqs = fastq_1 ? (fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ]) : []
            // A BAM file does not say whether its reads were paired, so without FASTQ
            // files that comes from the 'single_end' column, which defaults to paired end
            def single_end = fastq_1 ? !fastq_2 : (meta.single_end instanceof Boolean ? meta.single_end : false)
            def sample_meta = [ id: meta.id, condition: meta.condition, strandedness: meta.strandedness, single_end: single_end ]
            return [ meta.id, sample_meta, fastqs, genome_bam ?: null, transcriptome_bam ?: null, salmon_results ?: null ]
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .set { ch_samplesheet }

    // The modules that read the sheet themselves only need the sample and condition, once
    // per sample, in samplesheet order, whatever input columns the sheet has
    ch_samplesheet_file = channel
        .fromList(samplesheet_rows.collect { row -> row[0] }.unique { meta -> meta.id })
        .collectFile(name: 'samplesheet.csv', seed: 'sample,condition', sort: 'index', newLine: true) { meta ->
            "${meta.id},${meta.condition}"
        }
        .first()

    //
    // Create channel from the contrasts file provided through params.contrasts
    //
    // Validated here next to the samplesheet, so each sheet is read once and both are
    // emitted the same way: the rows for the consumers that work from the values, the
    // validated path for the processes that read the sheet themselves.
    //
    // Both the id, variable, reference, target[, blocking] columns and the deprecated
    // contrast, treatment, control columns are accepted. Rows are normalised to the
    // latter, which the downstream code and the modules reading the sheet use.
    //

    // Fails early with a clear error when the contrasts file does not exist
    file(contrasts, checkIfExists: true)

    // Every column of the contrastsheet is a meta field, so a row comes back as a
    // single element list holding the meta map, not as the map itself
    def contrast_metas = samplesheetToList(contrasts, "${projectDir}/assets/schema_contrasts.json").collect { row -> row[0] }

    if (contrast_metas.any { meta -> isLegacyContrastRow(meta) }) {
        log.warn("The contrastsheet columns 'contrast', 'treatment' and 'control' are deprecated and will be removed in a future release. Please use 'id', 'target' and 'reference' instead, together with a 'variable' column set to 'condition'. The columns map as contrast -> id, treatment -> target and control -> reference")
    }
    if (contrast_metas.any { meta -> meta.blocking }) {
        log.warn("The contrastsheet 'blocking' column is ignored, blocking factors are not supported")
    }

    def contrast_rows = contrast_metas.collect { meta ->
        def contrast_row = normaliseContrastRow(meta)
        validateInputContrastsheet([[contrast_row]])
        return contrast_row
    }

    ch_contrastsheet = channel.fromList(contrast_rows)

    // The modules that read the sheet themselves expect the legacy header, so they get
    // a canonical copy of the normalised rows
    ch_contrastsheet_file = channel
        .fromList(contrast_rows)
        .collectFile(name: 'contrastsheet.csv', seed: 'contrast,treatment,control', sort: 'index', newLine: true) { row ->
            "${row.contrast},${row.treatment},${row.control}"
        }
        .first()

    emit:
    samplesheet        = ch_samplesheet                                      // channel: [ val(meta), [ fastq ], genome_bam, transcriptome_bam, salmon_results ], absent inputs null
    samplesheet_file   = ch_samplesheet_file                                 // channel: path(samplesheet.csv), with the sample, condition header, one row per sample
    contrastsheet      = ch_contrastsheet                                    // channel: [ contrast, treatment, control ]
    contrastsheet_file = ch_contrastsheet_file                               // channel: path(contrastsheet.csv), with the contrast, treatment, control header
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
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
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

    if (params.source) {
        log.warn("The '--source' parameter is deprecated and ignored: the input type is worked out from the samplesheet columns. It will be removed in a future release.")
    }
}

//
// Which input types the samplesheet gives, from the columns that are filled. Methods
// compare the samples of a contrast, so a column is filled for every sample or for none
//
def getSamplesheetInputTypes(samplesheet_rows) {
    def columns = [ fastq: [1, 'fastq_1'], genome_bam: [3, 'genome_bam'], transcriptome_bam: [4, 'transcriptome_bam'], salmon_results: [5, 'salmon_results'] ]
    def samples = samplesheet_rows.collect { row -> row[0].id }.unique()

    def input_types = columns.collectEntries { input_type, column ->
        def (index, name) = column
        def samples_with = samplesheet_rows.findAll { row -> row[index] }.collect { row -> row[0].id }.unique()
        def samples_without = samples - samples_with
        if (samples_with && samples_without) {
            error("Please check input samplesheet -> The '${name}' column must be filled for every sample or for none, as the samples of a contrast are compared with each other. It is empty for: ${samples_without.join(', ')}")
        }
        return [ (input_type): samples_with as boolean ]
    }

    if (!input_types.any { _input_type, given -> given }) {
        error("Please check input samplesheet -> Each sample needs at least one of the 'fastq_1', 'genome_bam', 'transcriptome_bam' or 'salmon_results' columns.")
    }
    if (input_types.transcriptome_bam && input_types.salmon_results) {
        error("Please check input samplesheet -> The 'transcriptome_bam' and 'salmon_results' columns cannot both be given, as either one is quantified for the transcript level analyses.")
    }

    return input_types
}

//
// Input types of the samplesheet at --input, for the callers that need them as values
// to decide which steps run
//
def samplesheetInputRoutes(input) {
    return getInputRoutes(getSamplesheetInputTypes(samplesheetToList(input, "${projectDir}/assets/schema_input.json")))
}

//
// How each part of the pipeline gets its input. A BAM file or Salmon result given in the
// samplesheet is used instead of making it from the FASTQ files, which are only processed
// for what nothing else in the samplesheet provides
//
def getInputRoutes(input_types) {
    def star = input_types.fastq && !input_types.genome_bam && !params.skip_alignment && params.aligner in ['star', 'star_salmon']
    def star_salmon = input_types.transcriptome_bam || (star && params.aligner == 'star_salmon')
    def pseudo_salmon = input_types.fastq && !input_types.transcriptome_bam && !input_types.salmon_results && params.pseudo_aligner == 'salmon'

    return [
        input_types:    input_types,
        fastq:          star || pseudo_salmon,                 // FASTQ files are trimmed and quality checked
        star:           star,                                  // FASTQ files are aligned with STAR
        genome:         input_types.genome_bam || star,        // genome BAM files reach the exon and junction methods
        star_salmon:    star_salmon,                           // transcriptome BAM files are quantified with Salmon
        pseudo_salmon:  pseudo_salmon,                         // FASTQ files are quantified with Salmon
        salmon_results: input_types.salmon_results as boolean, // Salmon results are taken from the samplesheet
    ]
}

//
// Check that the input routes can run the methods asked for
//
def validateInputRoutes(input_routes) {
    if (input_routes.input_types.fastq && !input_routes.fastq) {
        log.info("The FASTQ files of the samplesheet are not used: BAM files or Salmon results are given for every analysis that would use them.")
    }

    // ASpli reads genome BAM files, which only the STAR alignment or a genome BAM column provide
    if (params.aspli && !input_routes.genome) {
        error("--aspli needs genome BAM files: give a 'genome_bam' samplesheet column, or FASTQ files without '--skip_alignment'.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs, genome_bams, transcriptome_bams, salmon_results) = input[1..5]
    def sample = metas[0].id

    // Check that multiple runs of the same sample are of the same strandedness
    def strandedness_ok = metas.collect { meta -> meta.strandedness }.unique().size == 1
    if (!strandedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same strandedness: ${sample}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${sample}")
    }

    if (metas.collect { meta -> meta.condition }.unique().size() > 1) {
        error("Please check input samplesheet -> Multiple rows of a sample must have the same condition: ${sample}")
    }

    if (fastqs.any { run -> run } && !fastqs.every { run -> run }) {
        error("Please check input samplesheet -> Multiple rows of a sample must all have FASTQ files or none: ${sample}")
    }

    [ genome_bam: genome_bams, transcriptome_bam: transcriptome_bams, salmon_results: salmon_results ].each { column, values ->
        if (values.collect { value -> value?.toString() }.unique().size() > 1) {
            error("Please check input samplesheet -> Multiple rows of a sample must give the same '${column}': ${sample}")
        }
    }

    return [ metas[0], fastqs.flatten(), genome_bams[0], transcriptome_bams[0], salmon_results[0] ]
}

//
// Whether a contrastsheet row uses the deprecated contrast, treatment, control columns
//
def isLegacyContrastRow(meta) {
    return ['contrast', 'treatment', 'control'].any { key -> meta[key] }
}

//
// Normalise a contrastsheet row to [contrast, treatment, control]. The columns map as
// id -> contrast, target -> treatment and reference -> control, so the log fold change
// is target over reference. The schema has already checked that a row uses one complete
// set of columns, so the values are present and the sets are not mixed
//
def normaliseContrastRow(meta) {
    if (isLegacyContrastRow(meta)) {
        return [ contrast: meta.contrast, treatment: meta.treatment, control: meta.control ]
    }
    return [ contrast: meta.id, treatment: meta.target, control: meta.reference ]
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
// Exit pipeline if rMATS requested with mixed single and paired end samples
//
def rmatsReadError(samples) {
    samples
        .map { row -> row[0].single_end }
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
        .map { row -> row[0].strandedness }
        .unique()
        .collect()
        .map {
            if (it.size() > 1) {
                error("Please check input samplesheet -> Cannot run rMats with mixed stranded samples.")
            }
        }
}

//
// Fail when the reference sequences in a BAM header, by name and length, are not all in
// the index of the FASTA file the pipeline uses for them, as a BAM file aligned to another
// reference would be quantified against the wrong sequences. The stub header is a
// placeholder, so stub runs are not checked
//
def checkBamReferenceSequences(meta, sequences_file, fai_file, reference_name) {
    if (workflow.stubRun) {
        return
    }

    def reference_lengths = [:]
    fai_file.eachLine { line ->
        def fields = line.split('\t')
        reference_lengths[fields[0]] = fields[1]
    }

    def bam_sequences = sequences_file.readLines().findAll { line -> line.startsWith('@SQ') }.collect { line ->
        def tags = line.split('\t').drop(1).collectEntries { tag -> [ (tag.take(2)): tag.drop(3) ] }
        [ name: tags.SN, length: tags.LN ]
    }
    if (!bam_sequences) {
        error("The BAM file of sample '${meta.id}' has no @SQ lines in its header, so it cannot be checked against ${reference_name}.")
    }

    def mismatched = bam_sequences.findAll { sequence -> reference_lengths[sequence.name] != sequence.length }
    if (mismatched) {
        def examples = mismatched.take(5).collect { sequence -> "${sequence.name} (${sequence.length} bp)" }.join(', ')
        error(
            "The BAM file of sample '${meta.id}' was not aligned to ${reference_name}: ${mismatched.size()} of its ${bam_sequences.size()} reference sequences are missing from ${fai_file.name} or have another length, for example ${examples}. " +
            "Give the pipeline the references the BAM files were aligned to. For BAM files from nf-core/rnaseq that is its genome FASTA and, for transcriptome BAM files, the transcript FASTA it wrote to '<outdir>/genome', given with '--transcript_fasta'."
        )
    }
}

//
// Median read length of the reads of all the samples, from the read length histograms
// ('RL' lines) of their samtools stats files, given as [ meta, stats ] pairs. rMATS takes
// one read length per run and recommends the median when the reads vary in length.
// Returns the median, which is null when no sample has reads (as in a stub run), and the
// ids of the samples with no reads
//
def getSamtoolsStatsMedianReadLength(samples_stats) {
    def counts = [:].withDefault { 0L }
    def samples_without_reads = []
    samples_stats.each { meta, stats_file ->
        def sample_reads = 0L
        stats_file.eachLine { line ->
            def fields = line.split('\t')
            if (fields[0] == 'RL') {
                counts[fields[1].toInteger()] += fields[2].toLong()
                sample_reads += fields[2].toLong()
            }
        }
        if (sample_reads == 0L) {
            samples_without_reads << meta.id
        }
    }
    def total = counts.values().sum() ?: 0L
    def seen = 0L
    def median = counts.keySet().sort().find { length ->
        seen += counts[length]
        seen * 2 >= total
    }
    return [median: median, samples_without_reads: samples_without_reads.sort()]
}

//
// Exit pipeline if ASpli requested with mixed single and paired end samples
//
def aspliReadError(samples) {
    samples
        .map { row -> row[0].single_end }
        .unique()
        .collect()
        .map {
            if (it.size() > 1) {
                error("Please check input samplesheet -> Cannot run ASpli with mixed single and paired end samples.")
            }
        }
}

//
// Call a library forward, reverse or unstranded from the reads assigned to each strand.
// Copied from nf-core/rnaseq 3.22.2:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
//
def calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def totalFragments = forwardFragments + reverseFragments + unstrandedFragments
    def totalStrandedFragments = forwardFragments + reverseFragments

    def strandedness = 'undetermined'
    if (totalStrandedFragments > 0) {
        def forwardProportion = forwardFragments / (totalStrandedFragments as double)
        def reverseProportion = reverseFragments / (totalStrandedFragments as double)
        def proportionDifference = Math.abs(forwardProportion - reverseProportion)

        if (forwardProportion >= stranded_threshold) {
            strandedness = 'forward'
        }
        else if (reverseProportion >= stranded_threshold) {
            strandedness = 'reverse'
        }
        else if (proportionDifference <= unstranded_threshold) {
            strandedness = 'unstranded'
        }
    }

    // Without any fragment, as from a stub run, the percentages are 0 rather than NaN,
    // which cannot be written to JSON
    def scale = totalFragments > 0 ? 100 / (totalFragments as double) : 0

    return [
        inferred_strandedness: strandedness,
        forwardFragments: forwardFragments * scale,
        reverseFragments: reverseFragments * scale,
        unstrandedFragments: unstrandedFragments * scale,
    ]
}

//
// Parse the Salmon 'lib_format_counts.json' file to get the inferred strandedness.
// Copied from nf-core/rnaseq 3.22.2:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
//
def getSalmonInferredStrandedness(json_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    // The Salmon stub writes an empty file, which is read as no counts at all
    def json_text = json_file.text
    def libCounts = json_text.trim() ? new groovy.json.JsonSlurper().parseText(json_text) : [:]

    def forwardKeys = ['SF', 'ISF', 'MSF', 'OSF']
    def reverseKeys = ['SR', 'ISR', 'MSR', 'OSR']

    // Every fragment Salmon assigns has a strand, even in an unstranded library, so these
    // counts are in practice always zero: https://groups.google.com/g/sailfish-users/c/yxzBDv6NB6I
    def unstrandedKeys = ['IU', 'U', 'MU']

    def forwardFragments = forwardKeys.collect { key -> libCounts[key] ?: 0 }.sum()
    def reverseFragments = reverseKeys.collect { key -> libCounts[key] ?: 0 }.sum()
    def unstrandedFragments = unstrandedKeys.collect { key -> libCounts[key] ?: 0 }.sum()

    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

//
// Parse the RSeQC 'infer_experiment.py' report to get the inferred strandedness, for
// both the single end and the paired end patterns. Copied from nf-core/rnaseq 3.22.2:
// https://github.com/nf-core/rnaseq/blob/3.22.2/subworkflows/local/utils_nfcore_rnaseq_pipeline/main.nf
//
def getInferexperimentStrandedness(inferexperiment_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def forwardFragments = 0
    def reverseFragments = 0
    def unstrandedFragments = 0

    inferexperiment_file.eachLine { line ->
        def unstranded_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
        def se_sense_matcher = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
        def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
        def pe_sense_matcher = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
        def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/

        if (unstranded_matcher) unstrandedFragments = unstranded_matcher[0][1].toFloat() * 100
        if (se_sense_matcher) forwardFragments = se_sense_matcher[0][1].toFloat() * 100
        if (se_antisense_matcher) reverseFragments = se_antisense_matcher[0][1].toFloat() * 100
        if (pe_sense_matcher) forwardFragments = pe_sense_matcher[0][1].toFloat() * 100
        if (pe_antisense_matcher) reverseFragments = pe_antisense_matcher[0][1].toFloat() * 100
    }

    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

//
// Get the strandedness from the library type Salmon detected with '--libType A', as written
// to its 'meta_info.json'. In alignment mode Salmon writes no 'lib_format_counts.json', so
// the call is Salmon's own rather than one made with the thresholds above
//
def getSalmonLibraryTypeStrandedness(json_file) {
    def json_text = json_file.text
    def meta_info = json_text.trim() ? new groovy.json.JsonSlurper().parseText(json_text) : [:]
    def library_type = meta_info.library_types ? meta_info.library_types[0] : ''

    def strandedness = 'undetermined'
    if (library_type in ['SF', 'ISF', 'MSF', 'OSF']) {
        strandedness = 'forward'
    }
    else if (library_type in ['SR', 'ISR', 'MSR', 'OSR']) {
        strandedness = 'reverse'
    }
    else if (library_type in ['U', 'IU', 'MU', 'OU']) {
        strandedness = 'unstranded'
    }

    return [inferred_strandedness: strandedness, library_type: library_type]
}
