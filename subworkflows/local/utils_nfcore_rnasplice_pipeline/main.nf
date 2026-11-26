//
// Subworkflow with functionality specific to the nf-core/rnasplice pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsHelp                } from 'plugin/nf-schema'
include { paramsSummaryLog          } from 'plugin/nf-schema'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { validateParameters        } from 'plugin/nf-schema'
include { checkCondaChannels        } from 'plugin/nf-core-utils'
include { citationsOnTheFly         } from 'plugin/nf-core-utils'
include { toolCitationText          } from 'plugin/nf-core-utils'
include { toolBibliographyText      } from 'plugin/nf-core-utils'
include { checkConfigProvided       } from 'plugin/nf-core-utils'
include { checkProfileProvided      } from 'plugin/nf-core-utils'
include { completionEmail           } from 'plugin/nf-core-utils'
include { completionSummary         } from 'plugin/nf-core-utils'
include { dumpParametersToJSON      } from 'plugin/nf-core-utils'
include { getWorkflowVersion        } from 'plugin/nf-core-utils'

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
    input             //  string: Path to input samplesheet, for the blocking column check
    contrasts         //  string: Path to the contrasts sheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    if (version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion(workflow.manifest.version, workflow.commitId)}")
        System.exit(0)
    }
    if (outdir) {
        dumpParametersToJSON(outdir.toString(), params)
    }
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba'])) {
        checkCondaChannels()
    }

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

    if (help || help_full) {
        log.info paramsHelp(
            [beforeText: before_text, afterText: after_text, command: command, showHidden: show_hidden, fullHelp: help_full],
            (help instanceof String && help != 'true') ? help : ''
        )
        exit 0
    }

    log.info before_text
    log.info paramsSummaryLog(workflow)
    log.info after_text
    if (validate_params) {
        validateParameters()
    }

    //
    // Check config provided to the pipeline
    //
    checkConfigProvided()
    checkProfileProvided(nextflow_cli_args)

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from the contrasts file provided through params.contrasts
    //
    // Both the id, variable, reference, target[, blocking] columns and the deprecated
    // contrast, treatment, control columns are accepted. Rows are normalised to the
    // latter, which the downstream code and the modules reading the sheet use, plus
    // blocking when it is set.
    //

    // Fails early with a clear error when the contrasts file does not exist
    file(contrasts, checkIfExists: true)

    // Every column of the contrastsheet is a meta field, so a row comes back as a
    // single element list holding the meta map, not as the map itself
    def contrast_metas = samplesheetToList(contrasts, "${projectDir}/assets/schema_contrasts.json").collect { row -> row[0] }

    if (contrast_metas.any { meta -> isLegacyContrastRow(meta) }) {
        log.warn("The contrastsheet columns 'contrast', 'treatment' and 'control' are deprecated and will be removed in a future release. Please use 'id', 'target' and 'reference' instead, together with a 'variable' column set to 'condition'. The columns map as contrast -> id, treatment -> target and control -> reference")
    }

    def contrast_rows = contrast_metas.collect { meta ->
        def contrast_row = normaliseContrastRow(meta)
        validateInputContrastsheet([[contrast_row]])
        return contrast_row
    }

    if (contrast_rows.any { row -> row.blocking }) {
        log.warn("The contrastsheet 'blocking' column is used by DEXSeq DEU (--dexseq_exon), DEXSeq DTU (--dexseq_dtu), edgeR DEU (--edger_exon), LeafCutter (--leafcutter) and IsoformSwitchAnalyzeR (--isoformswitchanalyzer) only, the other methods ignore it")
    }

    // Fails before any analysis when a blocking column is not in the samplesheet, as the
    // DEXSeq modules would only find out after the alignment and counting. The header is
    // read from the file, as the parsed samplesheet keeps only the columns of its schema
    def blocked_contrasts = contrast_rows.findAll { row -> row.blocking }
    if (blocked_contrasts) {
        def samplesheet_columns = file(input).withReader { reader -> reader.readLine() }.tokenize(',').collect { column -> column.trim().replaceAll(/^"|"$/, '') }
        blocked_contrasts.each { row ->
            def missing_columns = row.blocking.tokenize(';') - samplesheet_columns
            if (missing_columns) {
                error("Contrast '${row.contrast}': blocking columns missing from the samplesheet: ${missing_columns.join(', ')}")
            }
        }
    }

    ch_contrastsheet = channel.fromList(contrast_rows)

    emit:
    contrastsheet      = ch_contrastsheet                                    // channel: [ contrast, treatment, control, blocking (when set) ]
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
    def columns = [ fastq: 'fastq_1', genome_bam: 'genome_bam', transcriptome_bam: 'transcriptome_bam', salmon_results: 'salmon_results' ]
    def samples = samplesheet_rows.collect { row -> row.meta.id }.unique()

    def input_types = columns.collectEntries { input_type, name ->
        def samples_with = samplesheet_rows.findAll { row -> row[name] }.collect { row -> row.meta.id }.unique()
        def samples_without = samples - samples_with
        if (samples_with && samples_without) {
            error("Please check input samplesheet -> The '${name}' column must be filled for every sample or for none, as the samples of a contrast are compared with each other. It is empty for: ${samples_without.join(', ')}")
        }
        return [ (input_type): samples_with as boolean ]
    }

    def fastq_samples_without_strandedness = samplesheet_rows.findAll { row -> row.fastq_1 && !row.meta.strandedness }.collect { row -> row.meta.id }.unique()
    if (fastq_samples_without_strandedness) {
        error("Please check input samplesheet -> The 'strandedness' column must be filled for samples with FASTQ files, use 'auto' to infer it. It is empty for: ${fastq_samples_without_strandedness.join(', ')}")
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
// The samples of the samplesheet and how each part of the pipeline gets its input. These
// are plain values, not channels, as the input routes are tested in workflow `if`
// conditions. A sample may have one row per sequencing run: its runs are merged under its
// id, keeping the FASTQ files of all of them for concatenation, while the BAM files and
// Salmon results, which describe the whole sample, must be the same on every row
//
def parseSamplesheet(input) {
    def samplesheet_rows = readSamplesheet(input)
    def input_routes = getInputRoutes(getSamplesheetInputTypes(samplesheet_rows))
    validateInputRoutes(input_routes)

    def samples = samplesheet_rows
        .collect { row ->
            def meta = row.meta
            def fastqs = row.fastq_1 ? [ row.fastq_1, row.fastq_2 ].findAll() : []
            // A BAM file does not say whether its reads were paired, so without FASTQ
            // files that comes from the 'single_end' column
            def single_end = row.fastq_1 ? !row.fastq_2 : meta.single_end
            // BAM files and Salmon results carry no strandedness to infer from their reads, so
            // without FASTQ files a missing value is taken as unstranded
            def strandedness = meta.strandedness ?: 'unstranded'
            def sample_meta = [ id: meta.id, condition: meta.condition, strandedness: strandedness, single_end: single_end ]
            return [ meta: sample_meta, fastqs: fastqs, genome_bam: row.genome_bam ?: null, transcriptome_bam: row.transcriptome_bam ?: null, salmon_results: row.salmon_results ?: null ]
        }
        .groupBy { run -> run.meta.id }
        .collect { _id, runs -> validateInputSamplesheet(runs) }

    return [ samples: samples, input_routes: input_routes ]
}

//
// The validated samplesheet rows, each as a map of the meta map and the other columns by
// name. samplesheetToList returns the meta map, then the columns without a meta field in
// the order of the schema, which is where their names are taken from
//
def readSamplesheet(input) {
    def schema = "${projectDir}/assets/schema_input.json"
    def path_columns = new groovy.json.JsonSlurper().parse(file(schema)).items.properties
        .findAll { _name, property -> !property.meta }
        .keySet()
        .toList()
    return samplesheetToList(input, schema).collect { row ->
        [ meta: row[0] ] + [ path_columns, row.drop(1) ].transpose().collectEntries()
    }
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

    if (input_routes.star_salmon && input_routes.pseudo_salmon) {
        log.warn("Both --aligner=star_salmon and --pseudo_aligner=salmon specified. Downstream analyses will be performed on both salmon output files.")
    }

    // ASpli reads genome BAM files, which only the STAR alignment or a genome BAM column provide
    if (params.aspli && !input_routes.genome) {
        error("--aspli needs genome BAM files: give a 'genome_bam' samplesheet column, or FASTQ files without '--skip_alignment'.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(runs) {
    def metas = runs.collect { run -> run.meta }
    def fastqs = runs.collect { run -> run.fastqs }
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

    [ 'genome_bam', 'transcriptome_bam', 'salmon_results' ].each { column ->
        if (runs.collect { run -> run[column]?.toString() }.unique().size() > 1) {
            error("Please check input samplesheet -> Multiple rows of a sample must give the same '${column}': ${sample}")
        }
    }

    def first_run = runs[0]
    return [ first_run.meta, fastqs.flatten(), first_run.genome_bam, first_run.transcriptome_bam, first_run.salmon_results ]
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
    def row = [ contrast: meta.id, treatment: meta.target, control: meta.reference ]
    // As in nf-core/differentialabundance, an NA blocking value means no blocking
    return meta.blocking && meta.blocking != 'NA' ? row + [ blocking: meta.blocking ] : row
}

//
// The blocking variables of a contrast that go into the model of a method that fits the
// design itself. `sample_rows` are the samplesheet rows of the contrast samples, each with
// the sample's `condition`. A variable with one value in the contrast is left out with a
// warning. The pipeline stops when a variable is confounded with the condition, or when the
// variables leave no residual degrees of freedom, which the methods would only report as a
// rank-deficient design
//
def resolveBlockingVariables(contrast_id, method, blocking, sample_rows) {
    def kept = []
    blocking.tokenize(';').each { column ->
        def values = sample_rows.collect { row -> row[column]?.toString()?.trim() }
        if (values.any { value -> !value }) {
            error("${method} contrast '${contrast_id}': blocking column '${column}' has no value for every sample")
        }
        if (values.unique().size() == 1) {
            log.warn("${method} contrast '${contrast_id}': blocking column '${column}' has a single value in the samples of the contrast, so it is left out of the model")
            return
        }
        def conditions_by_value = sample_rows.groupBy { row -> row[column].toString().trim() }.collectEntries { value, rows -> [value, rows*.condition.unique()] }
        if (conditions_by_value.every { _value, conditions -> conditions.size() == 1 }) {
            error("${method} contrast '${contrast_id}': blocking column '${column}' is confounded with the condition, as every value of it occurs in one condition only")
        }
        kept << column
    }
    def n_conditions = sample_rows*.condition.unique().size()
    def n_parameters = n_conditions + kept.sum(0) { column -> sample_rows.collect { row -> row[column].toString().trim() }.unique().size() - 1 }
    if (kept && sample_rows.size() <= n_parameters) {
        error("${method} contrast '${contrast_id}': the blocking columns ${kept.join(', ')} leave no residual degrees of freedom with ${sample_rows.size()} samples")
    }
    return kept
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
// Citations of the tools that ran, from the versions their processes report and the
// meta.yml files of the modules, given in a fixed order as the first meta.yml naming a tool
// and the first version reported for it are the ones cited. A version names a tool by its
// conda package, which can differ from the meta.yml name by a 'bioconductor-' or 'r-'
// prefix, case and punctuation, so names are matched on citationToolKey. Citations are put
// on one line, as a line break would end the YAML block of the MultiQC methods section.
// Tools that report no version before the citations are built, such as MultiQC, which
// reads them, are given as extra_tools
//
def toolCitations(version_entries, meta_yml_paths, extra_tools) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def tool_names = [:]
    meta_yml_paths.each { path ->
        try {
            def tools = yaml.load(file(path).text)?.tools
            (tools instanceof List ? tools : [])
                .findAll { tool -> tool instanceof Map }
                .collectMany { tool -> tool.keySet() as List }
                .each { name -> tool_names.putIfAbsent(citationToolKey(name), name) }
        }
        catch (Exception e) {
            log.warn("The tools of ${path} could not be read for the MultiQC methods section: ${e.message}")
        }
    }

    def versions = version_entries
        .collectMany { entry -> entry instanceof Path ? readVersionsYml(entry) : [ entry ] }
        .findAll { _process, _tool, version -> version?.toString()?.trim() }
        .collect { process, tool, version -> [ process.toString(), tool_names[citationToolKey(tool)] ?: tool.toString(), version.toString() ] }
        .sort { version -> version.join('\t') }

    return citationsOnTheFly(versions, meta_yml_paths, extra_tools).collectEntries { tool, citation ->
        [ (tool): citation + [ citation: oneLine(citation.citation), bibliography: oneLine(citation.bibliography) ] ]
    }
}

//
// A tool name compared without case, punctuation or a 'bioconductor-' or 'r-' prefix
//
def citationToolKey(name) {
    return name.toString().toLowerCase().replaceFirst(/^(bioconductor|r)-/, '').replaceAll(/[^a-z0-9]/, '')
}

//
// The [ process, tool, version ] entries of a versions.yml file. The values are read as
// written, as a YAML load would turn an unquoted version such as 1.10 into the number 1.1
//
def readVersionsYml(path) {
    def root = new org.yaml.snakeyaml.Yaml().compose(new java.io.StringReader(path.text))
    if (!(root instanceof org.yaml.snakeyaml.nodes.MappingNode)) {
        return []
    }
    return root.value.collectMany { process ->
        if (!(process.valueNode instanceof org.yaml.snakeyaml.nodes.MappingNode)) {
            return []
        }
        return process.valueNode.value
            .findAll { tool -> tool.valueNode instanceof org.yaml.snakeyaml.nodes.ScalarNode }
            .collect { tool -> [ process.keyNode.value, tool.keyNode.value, tool.valueNode.value ] }
    }
}

//
// Text with its whitespace, line breaks included, collapsed to single spaces
//
def oneLine(text) {
    return text?.replaceAll(/\s+/, ' ')?.replace(' )', ')')
}

//
// MultiQC methods section, citing the given tools
//
def methodsDescriptionText(mqc_methods_yaml, citations) {
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

    meta["tool_citations"] = toolCitationText(citations)
    meta["tool_bibliography"] = toolBibliographyText(citations)

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
        .map { values ->
            if (values.size() > 1) {
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
        .map { values ->
            if (values.size() > 1) {
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
        .map { values ->
            if (values.size() > 1) {
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
