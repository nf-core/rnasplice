#!/usr/bin/env Rscript
# Scripts adjusted from F1000 workflow
# Please see following for details:
# https://f1000research.com/articles/7-952

# Parse command arguments

argv <- commandArgs(trailingOnly = TRUE)

argc <- length(argv)

name <- argv[1]

feature <- argv[2]

gene <- argv[3]

analysis <- argv[4]


# Attach required packages

library(stageR)


# Define helper functions

stripVersion <- function(x) {

    substr(x, 1, 15)

}

read.DEXSeq <- function(gene, feature) {

    # Read gene-level results

    resultsGene <- read.delim(gene)

    resultsGene <- setNames(resultsGene$padj, resultsGene$groupID)

    # Create a vector of screening hypothesis p-values

    pScreen <- resultsGene

    names(pScreen) <- stripVersion(names(pScreen))

    # Read feature-level results

    resultsFeature <- read.delim(feature)

    # Create a matrix of confirmation hypothesis p-values

    pConfirmation <- matrix(resultsFeature$pvalue, ncol = 1)

    dimnames(pConfirmation) <- list(
        stripVersion(resultsFeature$featureID),
        "transcript"
    )

    # Create a tx2gene table

    tx2gene <- resultsFeature[, c("featureID", "groupID"), drop = FALSE]

    tx2gene <- apply(tx2gene, 2, stripVersion)

    tx2gene <- as.data.frame(tx2gene)

    # Return created objects

    list(
        resultsGene    = resultsGene,
        resultsFeature = resultsFeature,
        pScreen        = pScreen,
        pConfirmation  = pConfirmation,
        tx2gene        = tx2gene
    )

}

read.DRIMSeq <- function(gene, feature) {

    # Read gene-level results

    resultsGene <- read.delim(gene)

    resultsGene <- setNames(resultsGene$padj, resultsGene$groupID)

    # Create a vector of screening hypothesis p-values

    pScreen <- resultsGene$pvalue

    names(pScreen) <- stripVersion(resultsGene$gene_id)

    # Read feature-level results

    resultsFeature <- read.delim(feature)

    # Create a matrix of confirmation hypothesis p-values

    pConfirmation <- matrix(resultsFeature$pvalue, ncol = 1)

    rownames(pConfirmation) <- stripVersion(resultsFeature$feature_id)

    # Create a tx2gene table

    tx2gene <- resultsFeature[, c("feature_id", "gene_id"), drop = FALSE]

    tx2gene <- apply(tx2gene, 2, stripVersion)

    tx2gene <- as.data.frame(tx2gene)

    # Return created objects

    list(
        resultsGene    = resultsGene,
        resultsFeature = resultsFeature,
        pScreen        = pScreen,
        pConfirmation  = pConfirmation,
        tx2gene        = tx2gene
    )

}


# Read analysis outputs

if (analysis == "dexseq") {

    output <- read.DEXSeq(
        gene    = gene,
        feature = feature
    )

} else if (analysis == "drimseq") {

    output <- read.DRIMSeq(
        gene    = gene,
        feature = feature
    )

}


# Create stageRTx object

object <- stageRTx(
    pScreen         = output[["pScreen"]],
    pConfirmation   = output[["pConfirmation"]],
    tx2gene         = output[["tx2gene"]],
    pScreenAdjusted = FALSE
)


# Adjust p-values in a two-stage analysis

object <- stageWiseAdjustment(
    object  = object,
    method  = "dtu",
    alpha   = 0.05,
    allowNA = TRUE
)

# Retrieve the stage-wise adjusted p-values

pvalue <- getAdjustedPValues(
    object               = object,
    onlySignificantGenes = FALSE,
    order                = FALSE
)

# Save objects to disk

saveRDS(
    object = object,
    file   = paste0("stageRTx.", name, ".rds")
)

saveRDS(
    object = pvalue,
    file = paste0("getAdjustedPValues.", name, ".rds")
)


# Save results to disk

write.table(
    x         = pvalue,
    file      = paste0("getAdjustedPValues.", name, ".tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)


# Print session information

citation("stageR")

sessionInfo()
