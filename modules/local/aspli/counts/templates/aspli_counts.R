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

#' Estimate the read length of a BAM file
#'
#' @param bam Path of the BAM file
#' @param probs Quantile of the read lengths to return
#' @param chunk Number of reads to read at a time
#'
#' @return Quantile of the query widths of all the primary mapped reads. The reads are
#'   read in chunks and only their lengths are tabulated, so memory does not grow with
#'   the size of the BAM file

estimate_read_length <- function(bam, probs = 0.05, chunk = 1000000L) {

    flag <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)

    param <- Rsamtools::ScanBamParam(what = "qwidth", flag = flag)

    bam_file <- Rsamtools::BamFile(bam, yieldSize = chunk)

    open(bam_file)

    on.exit(close(bam_file))

    counts <- integer(0)

    repeat {

        widths <- Rsamtools::scanBam(bam_file, param = param)[[1]][["qwidth"]]

        if (length(widths) == 0) break

        tabulated <- tabulate(widths)

        length(counts) <- max(length(counts), length(tabulated))

        counts[is.na(counts)] <- 0L

        counts[seq_along(tabulated)] <- counts[seq_along(tabulated)] + tabulated

    }

    if (sum(counts) == 0) {

        stop(paste0("No mapped reads to estimate the read length from in ", bam), call. = FALSE)

    }

    as.integer(which(cumsum(counts) / sum(counts) >= probs)[1])

}

#' Write a table of ASpli results with its row names as the first column
#'
#' @param x Data frame returned by an ASpli accessor
#' @param file Path of the TSV file to write
#' @param id Name of the column holding the row names

write_table <- function(x, file, id) {

    x <- data.frame(id = rownames(x), x, check.names = FALSE, row.names = NULL)

    colnames(x)[1] <- id

    write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

}


# Set defaults and parse the options provided via `task.ext.args`

opt <- list(
    max_intron_size    = 50000L,
    min_anchor         = 10L,
    junction_threshold = 5L
)

args_opt <- parse_args('$task.ext.args')

for (ao in names(args_opt)) {

    if (!ao %in% names(opt)) {

        stop(paste("Invalid option:", ao))

    }

    opt[[ao]] <- as(args_opt[[ao]], class(opt[[ao]]))

}

for (ao in names(opt)) {

    if (is.na(opt[[ao]]) || opt[[ao]] < 0) {

        stop(paste0("Option --", ao, " must be a non-negative integer."), call. = FALSE)

    }

}


# Read the inputs staged by Nextflow

prefix <- '${task.ext.prefix ?: meta.id}'

targets_file <- '$targets'

gtf <- '$gtf'

read_length <- as.integer('$read_length')

cores <- as.integer('$task.cpus')


# Read the targets table, one row per sample

samples <- read.delim(targets_file, colClasses = "character", check.names = FALSE)

required_columns <- c("sample", "condition", "bam", "single_end", "strandedness")

missing_columns <- setdiff(required_columns, colnames(samples))

if (length(missing_columns) > 0) {

    stop(paste0("Targets table is missing the columns: ", paste(missing_columns, collapse = ", ")), call. = FALSE)

}

if (anyDuplicated(samples[["sample"]]) > 0) {

    stop("Each sample of the targets table must be given once.", call. = FALSE)

}

missing_bams <- samples[["bam"]][!file.exists(samples[["bam"]])]

if (length(missing_bams) > 0) {

    stop(paste0("Missing BAM files: ", paste(missing_bams, collapse = ", ")), call. = FALSE)

}

if (length(unique(samples[["condition"]])) < 2) {

    stop("ASpli needs samples from at least two conditions.", call. = FALSE)

}


# Library type. ASpli counts every sample with one library type, as single reads
# (SE) or as read pairs (PE)

single_end <- unique(tolower(samples[["single_end"]]))

if (length(single_end) > 1) {

    stop("ASpli cannot count single end and paired end samples together.", call. = FALSE)

}

lib_type <- if (single_end == "true") "SE" else "PE"


# Strand mode. ASpli only uses it to give a strand to read pairs, and counts reads
# over bins and junctions whatever their strand, so a forward library maps to 1,
# a reverse library to 2 and an unstranded one to 0, as in GenomicAlignments

strand_modes <- c(unstranded = 0L, forward = 1L, reverse = 2L)

unknown_strandedness <- setdiff(samples[["strandedness"]], names(strand_modes))

if (length(unknown_strandedness) > 0) {

    stop(
        paste0(
            "Strandedness must be one of ", paste(names(strand_modes), collapse = ", "),
            ", got: ", paste(unknown_strandedness, collapse = ", ")
        ),
        call. = FALSE
    )

}

strandedness <- unique(samples[["strandedness"]])

strand_mode <- if (length(strandedness) == 1) strand_modes[[strandedness]] else 0L

if (length(strandedness) > 1) {

    message("The samples have different strandedness, so read pairs are given no strand (strandMode 0).")

}


# Read length. ASpli only counts a read over an exon-intron boundary when it overlaps
# it by at least this many bases, so a read shorter than it is never counted there. When
# it is not given it is the 5th percentile of the read lengths of the BAM file with the
# shortest reads, so that trimmed reads of varying length are still counted

read_length_source <- "given"

if (is.na(read_length) || read_length <= 0) {

    read_length <- min(vapply(samples[["bam"]], estimate_read_length, integer(1)))

    read_length_source <- "estimated"

}

message(paste0("Library type: ", lib_type, ", strand mode: ", strand_mode, ", read length: ", read_length, " (", read_length_source, ")"))


# Attach required packages

suppressPackageStartupMessages(library(ASpli))


# Bin the genome

txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(gtf, format = "gtf"))

features <- binGenome(txdb, logTo = paste0(prefix, ".binGenome.log"), cores = cores)


# Count reads over genes, bins and junctions. The order of the conditions is the order
# they first appear in the targets table, and is kept in the counts objects for the
# contrasts to be matched against

targets <- data.frame(
    row.names        = samples[["sample"]],
    bam              = samples[["bam"]],
    condition        = samples[["condition"]],
    stringsAsFactors = FALSE
)

gbcounts <- gbCounts(
    features      = features,
    targets       = targets,
    minReadLength = read_length,
    maxISize      = opt[["max_intron_size"]],
    minAnchor     = opt[["min_anchor"]],
    libType       = lib_type,
    strandMode    = strand_mode
)

jcounts <- jCounts(
    counts        = gbcounts,
    features      = features,
    minReadLength = read_length,
    threshold     = opt[["junction_threshold"]],
    minAnchor     = opt[["min_anchor"]],
    libType       = lib_type,
    strandMode    = strand_mode
)


# Export the R objects for the differential usage of each contrast

saveRDS(features, file = paste0(prefix, ".features.rds"))

saveRDS(gbcounts, file = paste0(prefix, ".gbcounts.rds"))

saveRDS(jcounts, file = paste0(prefix, ".jcounts.rds"))


# Export the counts and the junction usage tables

write_table(countsg(gbcounts), paste0(prefix, ".gene_counts.tsv"), "gene")

write_table(countsb(gbcounts), paste0(prefix, ".bin_counts.tsv"), "bin")

write_table(countsj(gbcounts), paste0(prefix, ".junction_counts.tsv"), "junction")

write_table(countse1i(gbcounts), paste0(prefix, ".e1i_counts.tsv"), "bin")

write_table(countsie2(gbcounts), paste0(prefix, ".ie2_counts.tsv"), "bin")

write_table(junctionsPIR(jcounts), paste0(prefix, ".junctions.pir.tsv"), "junction")

write_table(junctionsPJU(jcounts), paste0(prefix, ".junctions.pju.tsv"), "junction")

write_table(irPIR(jcounts), paste0(prefix, ".ir.pir.tsv"), "junction")

write_table(altPSI(jcounts), paste0(prefix, ".alt.psi.tsv"), "junction")

write_table(esPSI(jcounts), paste0(prefix, ".es.psi.tsv"), "junction")


# Export the settings the samples were counted with

settings <- data.frame(
    setting = c("library_type", "strand_mode", "read_length", "read_length_source", "max_intron_size", "min_anchor", "junction_threshold", "conditions"),
    value   = c(lib_type, strand_mode, read_length, read_length_source, opt[["max_intron_size"]], opt[["min_anchor"]], opt[["junction_threshold"]], paste(condition.order(gbcounts), collapse = ","))
)

write.table(settings, file = paste0(prefix, ".settings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


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
