#!/usr/bin/env Rscript
# Scripts adjusted from F1000 workflow
# Please see following for details:
# https://f1000research.com/articles/7-952

# Parse command arguments

argv <- commandArgs(trailingOnly = TRUE)

argc <- length(argv)

feature <- argv[1]

gene <- argv[2]

analysis <- argv[3]

# Attach required packages

library(stageR)

# Define helper functions

stripVersion <- function(x) {

    substr(x, 1, 15)

}

read.DEXSeq <- function(gene, feature) {

    # Read gene-level results

    resultsGene <- readRDS(gene)

    # Create a vector of screening hypothesis p-values

    pScreen <- resultsGene

    names(pScreen) <- stripVersion(names(pScreen))

    # Read feature-level results

    resultsFeature <- readRDS(feature)

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

    resultsGene <- readRDS(gene)

    # Create a vector of screening hypothesis p-values

    pScreen <- resultsGene$pvalue

    names(pScreen) <- stripVersion(resultsGene$gene_id)

    # Read feature-level results

    resultsFeature <- readRDS(feature)

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

print(feature)

print(gene)

print(analysis)

if (analysis == "dexseq") {

    outputs <- mapply(
        FUN = read.DEXSeq,
        gene = gene,
        feature = feature,
        SIMPLIFY = FALSE
    )

} else if (analysis == "drimseq") {

    outputs <- mapply(
        FUN = read.DRIMSeq,
        gene = gene,
        feature = feature,
        SIMPLIFY = FALSE
    )

}

# Create stageRTx object

objects <- mapply(
    FUN = stageRTx,
    pScreen = lapply(outputs, "[[", "pScreen"),
    pConfirmation = lapply(outputs, "[[", "pConfirmation"),
    tx2gene = lapply(outputs, "[[", "tx2gene"),
    MoreArgs = list(pscreenAdjusted = FALSE),
    SIMPLIFY = FALSE
)

# Adjust p-values in a two-stage analysis

objects <- lapply(
    X = objects,
    FUN = stageWiseAdjustment,
    method = "dtu",
    alpha = 0.05,
    allowNA = TRUE
)

# Retrieve the stage-wise adjusted p-values

pvalues <- lapply(
    X = objects,
    FUN = getAdjustedPValues,
    order = FALSE,
    onlySignificantGenes = FALSE
)

# Save objects to disk

mapply(
    saveRDS,
    object = objects,
    file = paste0("stageRTx.", names, ".rds")
)

mapply(
    saveRDS,
    object = pvalues,
    file = paste0("getAdjustedPValues", names, ".rds")
)

# Save results to disk

mapply(
    write.table,
    x = pvalues,
    file = paste0("getAdjustedPValues", names, ".tsv"),
    MoreArgs = list(sep = "\t", quote = FALSE, row.names = FALSE)
)

# Print session information

citation("stageR")

sessionInfo()
