#!/usr/bin/env Rscript

library(DEXSeq)
library(BiocParallel)

args <- commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 3) {

    stop("Usage: run_dexseq.R <sample.data> <d.counts> <ncores> <denominator>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

sample.data <- args[1] # sample.data from DRIMSeq run_drimseq_filter.R
d.counts    <- args[2] # d.counts from DRIMSeq run_drimseq_filter.R
ncores      <- args[3] # MultiCoreParam ncores

if (length(args) == 4) {

    denominator <- args[4]  # denominator for lfc set by user

} else {

    denominator <- ""       # denominator for lfc as default "" meaning is set as first sample condition

}

######################################
######## Define Multicoreparams ######
######################################

BPPARAM <- BiocParallel::MulticoreParam(ncores, stop.on.error = TRUE)

######################################
######## Run DEXseq analysis #########
######################################

# Read DRIMSeq sample data and counts
sample.data <- read.table(sample.data, sep = "\t", header = TRUE)
d.counts <- read.table(d.counts, sep = "\t", header = TRUE)

# set colnames of sample.data
colnames(sample.data) <- c("sample", "condition")

# Take count data from same filtered DRIMSeq
count.data <- round(d.counts[,(3:ncol(d.counts))])

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
                            featureID = d.counts$feature_id,
                            groupID = d.counts$gene_id)

dxd <- DEXSeq::estimateSizeFactors(dxd)

dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE, BPPARAM = BPPARAM)

# Looks for condition specific difference in tx proportions
dxd <- DEXSeq::testForDEU(dxd, reducedModel = reducedModel, BPPARAM = BPPARAM)

# Define sample col used for lfc calculation - at current this is fixed to required "condition" column
fitExpToVar <- "condition"

# Get fold changes based on fitExpToVar col in colData and denominator defines baseline for lfc
dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = fitExpToVar, denominator = denominator, BPPARAM = BPPARAM)

# Get Results
dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)

# Get q vals
qval <- DEXSeq::perGeneQValue(dxr)

# Format q vals
dxr.g <- data.frame(gene = names(qval),qval)

# dxr tsv
keep <- colnames(dxr)[sapply(dxr, class) %in% c("numeric", "character")]
dxr.tsv <- dxr[ , colnames(dxr) %in% keep]

################################
######### Save outputs #########
################################

# dxd
saveRDS(dxd, "dxd.rds")

# results
saveRDS(dxr, "dxr.rds")
write.table(dxr.tsv , "dxr.tsv", sep="\t", quote=FALSE, row.names = TRUE)

# qvals
saveRDS(qval, "qval.rds")
write.table(dxr.g, "dxr.g.tsv", sep="\t", quote=FALSE, row.names = TRUE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DEXSeq")
sessionInfo()
