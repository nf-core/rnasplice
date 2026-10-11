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

#' Write the top differentially used genes of a contrast to a PDF file
#'
#' @param results Table returned by topSpliceDGE
#' @param file Path of the PDF file to write
#' @param lrt DGELRT object returned by diffSpliceDGE
#' @param n Number of genes to plot

write_plot_splice <- function(results, file, lrt, n = 10) {

    ids <- head(results[["Geneid"]], n = n)

    pdf(file = file, width = 11.7, height = 8.3)

    on.exit(dev.off())

    if (length(ids) == 0) {

        plot.new()

        title(main = "No genes to plot")

        return(invisible(NULL))

    }

    for (i in ids) { plotSpliceDGE(lrt, geneid = i) }

    invisible(NULL)

}


# Set defaults and parse the options provided via `task.ext.args`

opt <- list(
    n_edger_plot = 10L
)

args_opt <- parse_args('$task.ext.args')

for (ao in names(args_opt)) {

    if (!ao %in% names(opt)) {

        stop(paste("Invalid option:", ao))

    }

    opt[[ao]] <- as(args_opt[[ao]], class(opt[[ao]]))

}

if (is.na(opt[["n_edger_plot"]]) || opt[["n_edger_plot"]] < 0) {

    stop("Option --n_edger_plot must be a non-negative integer.", call. = FALSE)

}


# Read the inputs staged by Nextflow

featurecounts <- "featurecounts"

samplesheet <- '$samplesheet'

contrastsheet <- '$contrastsheet'


# Attach required packages

library(edgeR)


# Read samples table

samples <- read.csv(samplesheet, check.names = FALSE)

if (!all(c("sample", "condition") %in% colnames(samples))) {

    stop("Samplesheet must contain 'sample' and 'condition' column headers.", call. = FALSE)

}

samples <- samples[, c("sample", "condition"), drop = FALSE]

samples <- unique(samples)

if (anyDuplicated(samples[["sample"]]) > 0) {

    stop("Each sample of the samplesheet must be assigned to a single condition.", call. = FALSE)

}


# Read contrasts table

contrastsheet <- read.csv(contrastsheet, check.names = FALSE)

if (!all(c("contrast", "treatment", "control") %in% colnames(contrastsheet))) {

    stop("Contrastsheet must contain 'contrast', 'treatment' and 'control' column headers.", call. = FALSE)

}

contrastsheet <- contrastsheet[, c("contrast", "treatment", "control"), drop = FALSE]

unknown_conditions <- setdiff(
    c(contrastsheet[["treatment"]], contrastsheet[["control"]]),
    samples[["condition"]]
)

if (length(unknown_conditions) > 0) {

    stop(
        paste0(
            "Conditions of the contrastsheet missing from the samplesheet: ",
            paste(unknown_conditions, collapse = ", ")
        ),
        call. = FALSE
    )

}


# Read featureCounts files

files <- file.path(featurecounts, paste0(samples[["sample"]], ".featureCounts.tsv"))

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {

    stop(
        paste0("Missing featureCounts files: ", paste(missing_files, collapse = ", ")),
        call. = FALSE
    )

}

data <- lapply(files, read.delim, comment.char = "#")

data <- Reduce(merge, data)


# Extract counts matrix

counts <- data[, -(1:6), drop = FALSE]

counts <- as.matrix(counts)

colnames(counts) <- samples[["sample"]]


# Extract genes annotation

genes <- data[, 1:6, drop = FALSE]


# Create DGEList object

DGEList <- DGEList(
    counts  = counts,
    samples = samples,
    group   = samples[["condition"]],
    genes   = genes
)


# Normalization

keep <- filterByExpr(DGEList, group = DGEList\$samples\$group)

DGEList <- DGEList[keep, , keep.lib.sizes = FALSE]

DGEList <- calcNormFactors(DGEList)


# Create design matrix

group <- factor(DGEList\$samples\$group)

design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)


# Create contrasts matrix

contrasts <- makeContrasts(
    contrasts = paste(contrastsheet[["treatment"]], contrastsheet[["control"]], sep = "-"),
    levels    = levels(group)
)

colnames(contrasts) <- contrastsheet[["contrast"]]

contrasts <- contrasts[, colSums(contrasts != 0) > 0, drop = FALSE]

if (ncol(contrasts) == 0) {

    stop("No contrast of the contrastsheet compares two different conditions.", call. = FALSE)

}


# Estimate dispersions by empirical Bayes

DGEList <- estimateDisp(DGEList, design)


# Fit log-linear model to count data

DGEGLM <- glmQLFit(DGEList, design)


# Test for differential exon expression

DGELRT.exprs <- mapply(
    FUN      = glmQLFTest,
    contrast = asplit(contrasts, MARGIN = 2),
    MoreArgs = list(glmfit = DGEGLM),
    SIMPLIFY = FALSE
)

results.exprs <- lapply(DGELRT.exprs, topTags, n = Inf, sort.by = "none")


# Test for differential exon usage

DGELRT.usage <- mapply(
    FUN      = diffSpliceDGE,
    contrast = asplit(contrasts, MARGIN = 2),
    MoreArgs = list(glmfit = DGEGLM, geneid = "Geneid", exonid = "Start"),
    SIMPLIFY = FALSE
)

results.usage <- list(
    simes = lapply(DGELRT.usage, topSpliceDGE, test = "Simes", number = Inf),
    gene  = lapply(DGELRT.usage, topSpliceDGE, test = "gene", number = Inf),
    exon  = lapply(DGELRT.usage, topSpliceDGE, test = "exon", number = Inf)
)


# Save objects to disk

saveRDS(DGEList, file = "DGEList.rds")

saveRDS(DGEGLM,  file = "DGEGLM.rds")

saveRDS(DGELRT.exprs, file = "DGELRT.exprs.rds")

saveRDS(DGELRT.usage, file = "DGELRT.usage.rds")


# Save results to disk

mapply(
    write.csv,
    x = results.exprs,
    file = paste0("contrast_", colnames(contrasts), ".exprs.csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)

for (test in names(results.usage)) {

    mapply(
        write.csv,
        x = results.usage[[test]],
        file = paste0("contrast_", colnames(contrasts), ".usage.", test, ".csv"),
        MoreArgs = list(quote = FALSE, row.names = FALSE)
    )

}


# Save plots to disk

for (test in names(results.usage)) {

    mapply(
        write_plot_splice,
        results = results.usage[[test]],
        file = paste0("contrast_", colnames(contrasts), ".usage.", test, ".pdf"),
        lrt = DGELRT.usage,
        MoreArgs = list(n = opt[["n_edger_plot"]])
    )

}


# Save the software versions to disk

writeLines(
    c(
        '"${task.process}":',
        paste0("    r-base: ", paste(R.version[["major"]], R.version[["minor"]], sep = ".")),
        paste0("    bioconductor-edger: ", as.character(packageVersion("edgeR")))
    ),
    "versions.yml"
)


# Print session information

citation("edgeR")

sessionInfo()
