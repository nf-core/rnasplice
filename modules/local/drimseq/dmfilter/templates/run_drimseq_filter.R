#!/usr/bin/env Rscript

# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

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


# Set defaults and parse the options provided via `task.ext.args`

opt <- list(
    min_samps_gene_expr    = 4,
    min_samps_feature_expr = 2,
    min_samps_feature_prop = 2,
    min_feature_expr       = 10,
    min_feature_prop       = 0.1,
    min_gene_expr          = 10
)

args_opt <- parse_args('$task.ext.args')

for (ao in names(args_opt)) {

    if (!ao %in% names(opt)) {

        stop(paste("Invalid option:", ao))

    }

    opt[[ao]] <- as(args_opt[[ao]], class(opt[[ao]]))

}


# Read the inputs staged by Nextflow

counts_file <- '$counts'

samplesheet <- '$samplesheet'


# Attach required packages

library(DRIMSeq)


# Read in sample sheet

samps <- read.csv(samplesheet, sep = ",", header = TRUE, check.names = FALSE)

if (!"sample" %in% colnames(samps) || !"condition" %in% colnames(samps)) {

    stop("Samplesheet must contain 'sample' and 'condition' column headers.", call. = FALSE)

}

# Take only sample and condition columns

samps <- samps[, c("sample", "condition")]

# Filter for unique rows based on sample name

samps <- samps[!duplicated(samps[, "sample"]), ]

# Change name of cols for DRIMSeq

colnames(samps) <- c("sample_id", "condition")


# Read the tximport transcript counts (scaledTPM or dtuScaledTPM): transcript ID, gene ID
# and one column per sample

cts <- read.csv(counts_file, sep = "\\t", header = TRUE, check.names = FALSE)

counts <- data.frame(
    gene_id    = cts[[2]],
    feature_id = cts[[1]],
    cts[, -(1:2), drop = FALSE],
    check.names = FALSE
)

# tximport keeps transcripts missing from the tx2gene as genes of their own, which are
# not part of the annotation, so they are left out

counts <- counts[counts\$gene_id != counts\$feature_id, ]

# Filter for txs with > 0 counts across all samples

counts <- counts[rowSums(counts[, (3:ncol(counts))]) > 0, ]


# Create DRIMSeq data set

d <- dmDSdata(counts = counts, samples = samps)

# Filter out genes and features with low expression

d <- dmFilter(
    x                      = d,
    min_samps_gene_expr    = opt[["min_samps_gene_expr"]],
    min_samps_feature_expr = opt[["min_samps_feature_expr"]],
    min_samps_feature_prop = opt[["min_samps_feature_prop"]],
    min_feature_expr       = opt[["min_feature_expr"]],
    min_feature_prop       = opt[["min_feature_prop"]],
    min_gene_expr          = opt[["min_gene_expr"]]
)

# Take pre-filtered sample data from DRIMSeq object

sample.data <- DRIMSeq::samples(d)

# Take count data

d.counts <- counts(d)


# Save object to disk

saveRDS(
    object = d,
    file   = "dmDSdata.rds"
)


# Save results to disk

write.table(
    x         = sample.data,
    file      = "samples.tsv",
    sep       = "\\t",
    quote     = FALSE,
    row.names = FALSE
)

write.table(
    x         = d.counts,
    file      = "counts.tsv",
    sep       = "\\t",
    quote     = FALSE,
    row.names = FALSE
)


# Save the software versions to disk

writeLines(
    c(
        '"${task.process}":',
        paste0("    r-base: ", paste(R.version[["major"]], R.version[["minor"]], sep = ".")),
        paste0("    bioconductor-drimseq: ", as.character(packageVersion("DRIMSeq")))
    ),
    "versions.yml"
)


# Print session information

citation("DRIMSeq")

sessionInfo()
