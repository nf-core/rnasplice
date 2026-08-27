#!/usr/bin/env Rscript

# Scripts adjusted from F1000 workflow
# Please see following for details:
# https://f1000research.com/articles/7-952

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

stripVersion <- function(x) {

    sub("\\\\..*", "", x)

}

read.DEXSeq <- function(gene, feature) {

    # Read gene-level results

    resultsGene <- read.delim(gene)

    resultsGene <- setNames(resultsGene\$padj, resultsGene\$groupID)

    # Create a vector of screening hypothesis p-values

    pScreen <- resultsGene

    names(pScreen) <- stripVersion(names(pScreen))

    # Read feature-level results

    resultsFeature <- read.delim(feature)

    # Create a matrix of confirmation hypothesis p-values

    pConfirmation <- matrix(resultsFeature\$pvalue, ncol = 1)

    dimnames(pConfirmation) <- list(
        stripVersion(resultsFeature\$featureID),
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

    resultsGene <- setNames(resultsGene\$padj, resultsGene\$groupID)

    # Create a vector of screening hypothesis p-values

    pScreen <- resultsGene\$pvalue

    names(pScreen) <- stripVersion(resultsGene\$gene_id)

    # Read feature-level results

    resultsFeature <- read.delim(feature)

    # Create a matrix of confirmation hypothesis p-values

    pConfirmation <- matrix(resultsFeature\$pvalue, ncol = 1)

    rownames(pConfirmation) <- stripVersion(resultsFeature\$feature_id)

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


# Set defaults and parse the options provided via `task.ext.args`

opt <- list(
    alpha                  = 0.05,
    method                 = "dtu",
    allow_na               = TRUE,
    only_significant_genes = FALSE,
    order                  = FALSE
)

args_opt <- parse_args('$task.ext.args')

for (ao in names(args_opt)) {

    if (!ao %in% names(opt)) {

        stop(paste("Invalid option:", ao))

    }

    opt[[ao]] <- as(args_opt[[ao]], class(opt[[ao]]))

}


# Read the inputs staged by Nextflow

if ("$task.ext.prefix" != "null") {

    name <- "$task.ext.prefix"

} else {

    name <- "$meta.id"

}

feature <- '$feature_tsv'

gene <- '$gene_tsv'

analysis <- '$analysis_type'


# Attach required packages

library(stageR)


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

} else {

    stop(paste("Invalid analysis type:", analysis))

}


# Create stageRTx object

object <- stageRTx(
    pScreen         = output[["pScreen"]],
    pConfirmation   = output[["pConfirmation"]],
    tx2gene         = output[["tx2gene"]],
    pScreenAdjusted = switch(analysis, "dexseq" = TRUE, "drimseq" = FALSE)
)


# Adjust p-values in a two-stage analysis

object <- stageWiseAdjustment(
    object  = object,
    method  = opt[["method"]],
    alpha   = opt[["alpha"]],
    allowNA = opt[["allow_na"]]
)

# Retrieve the stage-wise adjusted p-values

pvalue <- getAdjustedPValues(
    object               = object,
    onlySignificantGenes = opt[["only_significant_genes"]],
    order                = opt[["order"]]
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
    sep       = "\\t",
    quote     = FALSE,
    row.names = FALSE
)


# Save the software versions to disk

writeLines(
    c(
        '"${task.process}":',
        paste0("    stageR: ", as.character(packageVersion("stageR")))
    ),
    "versions.yml"
)


# Print session information

citation("stageR")

sessionInfo()
