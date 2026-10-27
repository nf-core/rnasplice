#!/usr/bin/env Rscript

# NOTE: This file is a Nextflow template. Nextflow interpolates the process
#       variables and treats the backslash as an escape character, therefore
#       R's list/data.frame operator must be written as a backslash followed
#       by a dollar sign, and every literal backslash must be doubled.


# Define helper functions

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x) {

    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]

    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones

    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })

    parsed_args <- structure(
        lapply(args_vals, function(x) x[2]),
        names = lapply(args_vals, function(x) x[1])
    )

    parsed_args[!is.na(parsed_args)]

}

#' Write a table of ASpli results with its row names as the first column
#'
#' @param x Data frame returned by an ASpli accessor
#' @param file Path of the TSV file to write
#' @param id Name of the column holding the row names, or NULL to leave them out

write_table <- function(x, file, id = NULL) {

    if (!is.null(id)) {

        x <- data.frame(id = rownames(x), x, check.names = FALSE, row.names = NULL)

        colnames(x)[1] <- id

    }

    write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

}


# Set defaults and parse the options provided via `task.ext.args`. The defaults are
# the ones of ASpli

opt <- list(
    min_gen_reads  = 10,
    min_bin_reads  = 5,
    min_rds        = 0.05,
    min_avg_counts = 5,
    strong_filter  = TRUE,
    bin_fc         = 3,
    bin_fdr        = 0.05,
    bin_inclusion  = 0.2,
    bjs_inclusion  = 10.3,
    bjs_fdr        = 0.01,
    a_inclusion    = 0.3,
    a_fdr          = 0.01,
    l_inclusion    = 0.3,
    l_fdr          = 0.01
)

args_opt <- parse_args('$task.ext.args')

for (ao in names(args_opt)) {

    if (!ao %in% names(opt)) {

        stop(paste("Invalid option:", ao))

    }

    opt[[ao]] <- as(args_opt[[ao]], class(opt[[ao]]))

}

for (ao in names(opt)) {

    if (is.na(opt[[ao]])) {

        stop(paste0("Option --", ao, " has an invalid value."), call. = FALSE)

    }

}


# Read the inputs staged by Nextflow

prefix <- '${task.ext.prefix ?: meta.id}'

contrast_name <- '${meta.id}'

treatment <- '${meta.treatment}'

control <- '${meta.control}'


# Attach required packages

suppressPackageStartupMessages(library(ASpli))


# Read the objects of ASPLI_COUNTS

features <- readRDS('$features')

gbcounts <- readRDS('$gbcounts')

jcounts <- readRDS('$jcounts')


# Build the contrast vector. ASpli weights each condition by the value at its position
# in the condition order of the counts, so the treatment is set to 1, the control to
# -1 and every other condition to 0. The log fold changes are treatment over control

conditions <- condition.order(gbcounts)

missing_conditions <- setdiff(c(treatment, control), conditions)

if (length(missing_conditions) > 0) {

    stop(
        paste0(
            "Contrast '", contrast_name, "': conditions with no sample in the counts: ",
            paste(missing_conditions, collapse = ", ")
        ),
        call. = FALSE
    )

}

contrast <- stats::setNames(rep(0, length(conditions)), conditions)

contrast[[treatment]] <- 1

contrast[[control]] <- -1


# Differential gene expression and bin usage

gb_du <- gbDUreport(
    gbcounts,
    minGenReads = opt[["min_gen_reads"]],
    minBinReads = opt[["min_bin_reads"]],
    minRds      = opt[["min_rds"]],
    contrast    = contrast
)


# Differential junction usage. With the strong filter ASpli drops every event with a
# junction that does not pass the filter, which can leave too few junctions to test.
# With more than 24 conditions ASpli estimates the dispersion from a random sample of
# the conditions outside the contrast, so the seed is set for reproducible results

set.seed(1)

j_du <- tryCatch(
    jDUreport(
        jcounts,
        minAvgCounts = opt[["min_avg_counts"]],
        contrast     = contrast,
        strongFilter = opt[["strong_filter"]]
    ),
    error = function(e) {
        stop(
            paste0(
                "ASpli jDUreport failed for contrast '", contrast_name, "': ", conditionMessage(e),
                if (opt[["strong_filter"]]) " Setting '--strong_filter false' in the ext.args of ASPLI_DU may help." else ""
            ),
            call. = FALSE
        )
    }
)


# Integrate the bin and junction signals

splicing_report <- splicingReport(gb_du, j_du, gbcounts)

integrated_signals <- integrateSignals(
    splicing_report,
    jcounts,
    bin.FC         = opt[["bin_fc"]],
    bin.fdr        = opt[["bin_fdr"]],
    bin.inclussion = opt[["bin_inclusion"]],
    bjs.inclussion = opt[["bjs_inclusion"]],
    bjs.fdr        = opt[["bjs_fdr"]],
    a.inclussion   = opt[["a_inclusion"]],
    a.fdr          = opt[["a_fdr"]],
    l.inclussion   = opt[["l_inclusion"]],
    l.fdr          = opt[["l_fdr"]]
)


# Export the tables

write_table(genesDE(gb_du), paste0(prefix, ".genes_de.tsv"), "gene")

write_table(binsDU(gb_du), paste0(prefix, ".bins_du.tsv"), "bin")

write_table(localej(j_du), paste0(prefix, ".locale_junctions.tsv"), "junction")

write_table(localec(j_du), paste0(prefix, ".locale_clusters.tsv"), "cluster")

write_table(anchorj(j_du), paste0(prefix, ".anchor_junctions.tsv"), "junction")

write_table(anchorc(j_du), paste0(prefix, ".anchor_clusters.tsv"), "junction")

write_table(jir(j_du), paste0(prefix, ".ir_junctions.tsv"), "junction")

write_table(jes(j_du), paste0(prefix, ".es_junctions.tsv"), "junction")

write_table(jalt(j_du), paste0(prefix, ".alt_junctions.tsv"), "junction")

write_table(binbased(splicing_report), paste0(prefix, ".splicing_report_bins.tsv"))

write_table(localebased(splicing_report), paste0(prefix, ".splicing_report_locale.tsv"))

write_table(anchorbased(splicing_report), paste0(prefix, ".splicing_report_anchor.tsv"))

write_table(signals(integrated_signals), paste0(prefix, ".integrated_signals.tsv"))


# Export the R objects

saveRDS(
    list(gb_du = gb_du, j_du = j_du, splicing_report = splicing_report, integrated_signals = integrated_signals),
    file = paste0(prefix, ".aspli_du.rds")
)


# Export the integrated signals as an HTML report. The event graphs would need BAM
# files merged by condition, so they are not drawn and the BAM files are never read,
# but the function still expects a table of them

merged_bams <- data.frame(bam = paste0(conditions, ".bam"), condition = conditions)

exportIntegratedSignals(
    integrated_signals,
    output.dir = "integrated_signals",
    sr         = splicing_report,
    counts     = gbcounts,
    features   = features,
    asd        = jcounts,
    mergedBams = merged_bams,
    makeGraphs = FALSE
)

# ASpli writes the report to a directory named after the contrasted conditions
report_dir <- list.dirs("integrated_signals", recursive = FALSE)

if (length(report_dir) != 1) {

    stop("ASpli did not write the integrated signals report.", call. = FALSE)

}

file.rename(report_dir, paste0(prefix, "_integrated_signals"))


# Save the software versions to disk

writeLines(
    c(
        '"${task.process}":',
        paste0("    r-base: ", paste(R.version[["major"]], R.version[["minor"]], sep = ".")),
        paste0("    bioconductor-aspli: ", as.character(packageVersion("ASpli")))
    ),
    "versions.yml"
)


# Print session information

sessionInfo()
