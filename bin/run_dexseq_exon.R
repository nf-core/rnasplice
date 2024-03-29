#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

# Read command arguments

argv <- commandArgs(trailingOnly = TRUE)

argc <- length(argv)


# Parse command arguments

counts <- argv[1]

annotation <- argv[2]

samples <- argv[3]

contrasts <- argv[4]

ntop <- as.integer(argv[5])


# Load required packages

library(BiocParallel)

library(DEXSeq)


# Define helper functions

splitByContrast <- function(object, contrast) {

    keep <- colData(object)$condition %in% contrast

    object <- object[, keep, drop = FALSE]

    colData(object)$condition <- factor(
        x = colData(object)$condition,
        levels = contrast
    )

    object

}

DEXSeq <- function(object) {

    object <- estimateSizeFactors(object)

    object <- estimateDispersions(object)

    object <- testForDEU(object)

    object <- estimateExonFoldChanges(object)

    object

}

keepValidColumns <- function(x) {

    keep <- colnames(x)[sapply(x, class) %in% c("numeric", "character")]

    x <- x[ , colnames(x) %in% keep, drop = FALSE]

    x

}

vectorToDataFrame <- function(x) {

    x <- data.frame(groupID = names(x), padj = unname(x))

    x

}

write.plotDEXSeq <- function(x, file, ntop = 10) {

    ind <- which(x$pvalue <= sort(x$pvalue)[ntop], arr.ind = TRUE)

    ids <- x$groupID[ind]

    pdf(file = file, width = 11.7, height = 8.3)

    for (i in ids) {

        plotDEXSeq(x, geneID = i)

    }

    dev.off()

}

# Read samples table

samples <- read.csv(samples, stringsAsFactors = TRUE, check.names = FALSE)

samples <- samples[, c("sample", "condition"), drop = FALSE]

samples <- unique(samples)

# Read contrasts table

contrasts <- read.csv(contrasts, check.names = FALSE)

contrasts <- contrasts[, c("contrast", "treatment", "control"), drop = FALSE]

contrasts <- unique(contrasts)


# List counts files

files <- paste0(counts, "/", samples$sample, ".clean.count.txt")


# Create DEXSeqDataSet object

object <- DEXSeqDataSetFromHTSeq(
    countfiles = files,
    sampleData = samples,
    flattenedfile = annotation
)


# Define contrasts

groups <- levels(colData(object)$condition)

names <- contrasts$contrast

contrasts <- data.frame(A = contrasts$treatment, B = contrasts$control)

contrasts <- contrasts[contrasts$A != contrasts$B, , drop = FALSE]

contrasts <- asplit(contrasts, MARGIN = 1)


# Create objects list

objects <- mapply(
    FUN = splitByContrast,
    contrast = contrasts,
    MoreArgs = list(object = object),
    SIMPLIFY = FALSE
)


# Run DEXSeq workflow

objects <- lapply(objects, DEXSeq)

names(objects) <- names


# Run DEXSeqResults

results.exon <- lapply(objects, DEXSeqResults)

names(results.exon) <- names


# Run perGeneQValue

results.gene <- lapply(results.exon, perGeneQValue)

names(results.gene) <- names


# Save objects to disk

mapply(
    saveRDS,
    object = objects,
    file = paste0("DEXSeqDataSet.", names, ".rds")
)

mapply(
    saveRDS,
    object = results.exon,
    file = paste0("DEXSeqResults.", names, ".rds")
)

mapply(
    saveRDS,
    object = results.gene,
    file = paste0("perGeneQValue.", names, ".rds")
)


# Save results to disk

mapply(
    write.csv,
    x = lapply(results.exon, keepValidColumns),
    file = paste0("DEXSeqResults.", names, ".csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
    write.csv,
    x = lapply(results.gene, vectorToDataFrame),
    file = paste0("perGeneQValue.", names, ".csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)

# Save plots to disk

mapply(
    write.plotDEXSeq,
    x = results.exon,
    file = paste0("plotDEXSeq.", names, ".pdf"),
    MoreArgs = list(ntop = ntop)
)

# Print session information

citation("DEXSeq")

sessionInfo()
