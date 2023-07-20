//
// This file holds several functions specific to the workflow/rnasplice.nf in the nf-core/rnasplice pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowRnasplice {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        genomeExistsError(params, log)


        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }

        if (!params.gtf && !params.gff) {
            Nextflow.error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
        }

        if (params.gtf) {
            if (params.gff) {
                gtfGffWarn(log)
            }
            if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
                ncbiGenomeWarn(log)
            }
            if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
                ucscGenomeWarn(log)
            }
        }

        if (params.transcript_fasta) {
            transcriptsFastaWarn(log)
        }

        if (!params.skip_alignment) {
            if (!valid_params['aligners'].contains(params.aligner)) {
                Nextflow.error("Invalid option: '${params.aligner}'. Valid options for '--aligner': '${valid_params['aligners'].join(', ')}'.")
            }
        } else {
            if (!params.pseudo_aligner) {
                Nextflow.error("--skip_alignment specified without --pseudo_aligner...please specify e.g. --pseudo_aligner '${valid_params['pseudoaligners'][0]}'.")
            }
            skipAlignmentWarn(log)
        }

        if (params.pseudo_aligner) {
            if (!valid_params['pseudoaligners'].contains(params.pseudo_aligner)) {
                Nextflow.error("Invalid option: '${params.pseudo_aligner}'. Valid options for '--pseudo_aligner': '${valid_params['pseudoaligners'].join(', ')}'.")
            } else {
                if (!(params.salmon_index || params.transcript_fasta || (params.fasta && (params.gtf || params.gff)))) {
                    Nextflow.error("To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff.")
                }
            }
        }

        //
        //  SUPPA parameter checks
        //

        if (params.clusterevents_local_event && !params.diffsplice_local_event) {
            Nextflow.error("--clusterevents_local_event specified without --diffsplice_local_event... please specify e.g. --diffsplice_local_event=true")
        }

        if (params.clusterevents_isoform && !params.diffsplice_isoform) {
            Nextflow.error("--clusterevents_isoform specified without diffsplice_isoform... please specify e.g. --diffsplice_isoform=true")
        }

        //
        // Warn of duplicated analyses from aligner and pseudo_aligner
        //

        if (!params.skip_alignment) {
            if (params.aligner == "star_salmon"  && params.pseudo_aligner == "salmon") {
                log.warn "Both --aligner=star_salmon and --pseudo_aligner=salmon specified. Downstream analyses will be performed on both salmon output files."
            }
        }

    }

    //
    // Function to check whether biotype field exists in GTF file
    //
    public static Boolean biotypeInGtf(gtf_file, biotype, log) {
        def hits = 0
        gtf_file.eachLine { line ->
            def attributes = line.split('\t')[-1].split()
            if (attributes.contains(biotype)) {
                hits += 1
            }
        }
        if (hits) {
            return true
        } else {
            log.warn "=============================================================================\n" +
                "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
                "  Biotype QC will be skipped to circumvent the issue below:\n" +
                "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
                "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
                "==================================================================================="
            return false
        }
    }

    //
    // Function to generate an error if contigs in genome fasta file > 512 Mbp
    //
    public static void checkMaxContigSize(fai_file, log) {
        def max_size = 512000000
        fai_file.eachLine { line ->
            def lspl  = line.split('\t')
            def chrom = lspl[0]
            def size  = lspl[1]
            if (size.toInteger() > max_size) {
                Nextflow.error(
                    "=============================================================================\n" +
                    "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
                    "  ${chrom}: ${size}\n\n" +
                    "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
                    "  Please see:\n" +
                    "  https://github.com/nf-core/rnaseq/issues/744\n" +
                    "============================================================================="
                )
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Create MultiQC tsv custom content from a list of values
    //
    public static String multiqcTsvFromList(tsv_data, header) {
        def tsv_string = ""
        if (tsv_data.size() > 0) {
            tsv_string += "${header.join('\t')}\n"
            tsv_string += tsv_data.join('\n')
        }
        return tsv_string
    }

    //
    // Generate methods description for MultiQC
    //
    public static String toolCitationText(params) {

        // TODO Optionally add in-text citation tools to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {

        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Print a warning if using GRCh38 assembly from igenomes.config
    //
    private static void ncbiGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "==================================================================================="
    }

    //
    // Print a warning if using a UCSC assembly from igenomes.config
    //
    private static void ucscGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "==================================================================================="
    }

    //
    // Print a warning if both GTF and GFF have been provided
    //
    private static void gtfGffWarn(log) {
        log.warn "=============================================================================\n" +
            "  Both '--gtf' and '--gff' parameters have been provided.\n" +
            "  Using GTF file as priority.\n" +
            "==================================================================================="
    }

    //
    // Print a warning if using '--transcript_fasta'
    //
    private static void transcriptsFastaWarn(log) {
        log.warn "=============================================================================\n" +
            "  '--transcript_fasta' parameter has been provided.\n" +
            "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" +
            "  Please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/753\n" +
            "==================================================================================="
    }

    //
    // Print a warning if --skip_alignment has been provided
    //
    private static void skipAlignmentWarn(log) {
        log.warn "=============================================================================\n" +
            "  '--skip_alignment' parameter has been provided.\n" +
            "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" +
            "==================================================================================="
    }

    //
    // Exit pipeline if rMATS requested with mixed single and paired end samples
    //
    public static void rmatsReadError(reads, log) {
        reads
            .map { meta, fastq -> meta.single_end }
            .unique()
            .collect()
            .map {
                if (it.size() > 1) {
                    Nextflow.error("Please check input samplesheet -> Cannot run rMats with mixed single and paired end samples.")
                }
            }
    }

    //
    // Exit pipeline if rMATS requested with mixed stranded samples
    //
    public static void rmatsStrandednessError(reads, log) {
        reads
            .map { meta, fastq -> meta.strandedness }
            .unique()
            .collect()
            .map {
                if (it.size() > 1) {
                    Nextflow.error("Please check input samplesheet -> Cannot run rMats with mixed stranded samples.")
                }
            }
    }

    //
    // Create variable to check if samples have one condition or multiple
    //

    public static Boolean isSingleCondition(samplesheet) {

        def reader = samplesheet.splitCsv(header: true)

        def conditions = []

        reader.each { row -> conditions << row.condition }

        def single_condition = conditions.unique().size() == 1

        return single_condition

    }

}
