#!/usr/bin/env Rscript

# Parse command arguments

argv <- commandArgs(trailingOnly=TRUE)

argc <- length(argv)

# Load required packages

library(DEXSeq)

file <- args[1]

n <- args[2]

# Plot read count data, expression and exon usage

object <- readRDS(file)

top_dxr <- unique(object$groupID[order(object$padj)])[1:n]

pdf(file = "dexseq_plot.pdf", width = 11.7, height = 8.3)

for (i in 1:n) {

    plotDEXSeq(
        object = object,
        geneID = top_dxr[i],
        legend = TRUE,
        displayTranscripts = TRUE,
        splicing = TRUE,
        cex.axis = 1,
        cex = 1,
        lwd = 2
    )

}

dev.off()

# Print session information

citation("DEXSeq")

sessionInfo()
