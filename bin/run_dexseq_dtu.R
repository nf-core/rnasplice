#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT


# Load required packages
library(BiocParallel)
library(parallel)
library(DEXSeq)

# Parse command arguments
argv <- commandArgs(trailingOnly = TRUE)
argc <- length(argv)
if (argc < 3) {
  stop("Usage: run_dexseq_dtu.R <samples_table> <contrasts_table> <counts_table> [ncores]")
}
samples <- argv[1]
contrasts <- argv[2]
counts <- argv[3]
ncores <- if (argc >= 4) as.integer(argv[4]) else 1

# Register parallel backend
if (ncores > 1) {
  BPPARAM <- MulticoreParam(workers = ncores)
} else {
  BPPARAM <- SerialParam()
}

#############################
## Define helper functions ##
#############################

splitByContrast <- function(object, contrast) {
  keep <- colData(object)$condition %in% contrast
  object <- object[, keep, drop = FALSE]
  colData(object)$condition <- factor(
    x = colData(object)$condition,
    levels = contrast
  )
  object
}

DEXSeqPipeline <- function(object, BPPARAM = BiocParallel::bpparam()) {
  object <- estimateSizeFactors(object)
  object <- estimateDispersions(object, BPPARAM = SerialParam())
  object <- testForDEU(object, BPPARAM = BPPARAM)
  object <- estimateExonFoldChanges(object, BPPARAM = BPPARAM)
  object
}

keepValidColumns <- function(x) {
  keep <- vapply(x, function(col) is.numeric(col) || is.character(col), logical(1))
  x <- x[, keep, drop = FALSE]
  x
}

vectorToDataFrame <- function(x) {
  x <- data.frame(groupID = names(x), padj = unname(x))
  x
}

########################
## Read samples table ##
########################

samples <- read.delim(samples, stringsAsFactors = TRUE, check.names = FALSE)
samples$sample <- as.character(samples$sample)
colnames(samples) <- c("sample", "condition")

##########################
## Read contrasts table ##
##########################

contrasts <- read.csv(contrasts, check.names = FALSE)
contrasts <- contrasts[, c("contrast", "treatment", "control"), drop = FALSE]

# Read counts table
counts <- read.table(counts, sep = "\t", header = TRUE, check.names = FALSE)

annotation <- data.frame(
  featureID = counts$feature_id,
  groupID = counts$gene_id
)

keep <- -c(1:2)
rownames(counts) <- counts$feature_id
counts <- counts[, keep, drop = FALSE]
counts <- round(counts)

# Reorder count matrix by sample table
stopifnot(all(samples$sample %in% colnames(counts)))
stopifnot(all(colnames(counts) %in% samples$sample))
counts <- counts[, as.character(samples$sample), drop = FALSE]
stopifnot(all(colnames(counts) == as.character(samples$sample)))

# Create DEXSeqDataSet object
dexseq_dataset <- DEXSeqDataSet(
  countData = counts,
  sampleData = samples,
  featureID = annotation$featureID,
  groupID = annotation$groupID
)

######################
## Define contrasts ##
######################

groups <- levels(colData(dexseq_dataset)$condition)
contrasts <- data.frame(A = contrasts$treatment, B = contrasts$control)
# We expect that the input contrasts table is already well
if (any(contrasts$A == contrasts$B)) {
  stop("Invalid contrasts found. Please check the contrasts table.")
}
# contrasts <- contrasts[contrasts$A != contrasts$B, , drop = FALSE]
contrasts <- asplit(contrasts, MARGIN = 1)
contrast_names <- sapply(contrasts, paste, collapse = "-")

# Create objects list
objects <- mapply(
  FUN = splitByContrast,
  contrast = contrasts,
  MoreArgs = list(object = dexseq_dataset),
  SIMPLIFY = FALSE
)

# Run DEXSeq workflow
objects <- lapply(objects, DEXSeqPipeline, BPPARAM = BPPARAM)
names(objects) <- contrast_names

# Run DEXSeqResults
results_exon <- lapply(objects, DEXSeqResults)
names(results_exon) <- contrast_names

# Run perGeneQValue
results_gene <- lapply(results_exon, perGeneQValue)
names(results_gene) <- contrast_names

##########################
## Save objects to disk ##
##########################

mapply(
  saveRDS,
  object = objects,
  file = paste0("DEXSeqDataSet.", contrast_names, ".rds")
)

mapply(
  saveRDS,
  object = results_exon,
  file = paste0("DEXSeqResults.", contrast_names, ".rds")
)

mapply(
  saveRDS,
  object = results_gene,
  file = paste0("perGeneQValue.", contrast_names, ".rds")
)


##########################
## Save results to disk ##
##########################

mapply(
  write.table,
  x = lapply(results_exon, keepValidColumns),
  file = paste0("DEXSeqResults.", contrast_names, ".tsv"),
  MoreArgs = list(quote = FALSE, sep = "\t", row.names = FALSE)
)

mapply(
  write.table,
  x = lapply(results_gene, vectorToDataFrame),
  file = paste0("perGeneQValue.", contrast_names, ".tsv"),
  MoreArgs = list(quote = FALSE, sep = "\t", row.names = FALSE)
)

###############################
## Print session information ##
###############################

citation("DEXSeq")
sessionInfo()
