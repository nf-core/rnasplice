#!/usr/bin/env Rscript

library(DEXSeq)
library(DRIMSeq)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 1) {
  
  stop("Usage: run_dexseq.R <drimseq_filter_rds>", call.=FALSE)
  
}

######################################
########### Collect inputs ###########
######################################

d <- args[1]  # d.rds object from run_drimseq_filter.R

######################################
######## Run DEXseq analysis #########
######################################

# Read r object
d <- readRDS(d)

# Take pre-filtered sample data from DRIMSeq object
sample.data <- DRIMSeq::samples(d)

# set colnames of sample.data
colnames(sample.data) <- c("sample", "condition")

# Take count data from same filtered DRIMSeq object
count.data <- round(DRIMSeq::counts(d)[,(3:ncol(DRIMSeq::counts(d)))])

# Define models
fullModel <- as.formula("~sample + exon + condition:exon")
reducedModel <- as.formula("~sample + exon")

# DEXseq made for exon level but modified for DTU
# design specifies "exon" but should be read as "transcript"
# See F1000 workflow for more details:
# https://f1000research.com/articles/7-952

dxd <- DEXSeq::DEXSeqDataSet(countData = count.data,
                             sampleData = sample.data,
                             design = fullModel,
                             featureID = counts(d)$feature_id,
                             groupID = counts(d)$gene_id)

dxd <- DEXSeq::estimateSizeFactors(dxd)

dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE)

# Looks for condition specific difference in tx proportions
dxd <- DEXSeq::testForDEU(dxd, reducedModel = reducedModel)

# Get Results
dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)

# Get q vals
qval <- DEXSeq::perGeneQValue(dxr)

# Format q vals
dxr.g <- data.frame(gene = names(qval),qval)

################################
######### Save outputs #########
################################

# dxd
saveRDS(dxd, "dxd.rds")

# results
saveRDS(dxr, "dxr.rds")
write.table(dxr, "dxr.tsv", sep="\t", quote=FALSE, row.names = TRUE)

# qvals
saveRDS(qval, "qval.rds")
write.table(dxr.g, "dxr.g.tsv", sep="\t", quote=FALSE, row.names = TRUE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DEXSeq")
sessionInfo()