#!/usr/bin/env Rscript

library(DEXSeq)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 5) {

  stop("Usage: run_dexseq_exon.R <countFiles_dir> <flattenedFile> <samplesheet> <read_method> <ncores> <denominator>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

countFiles_dir <- args[1] # count files
flattenedFile <- args[2]  # gff
samplesheet <- args[3]    # samplesheet
read_method <- args[4]    # either HTSeq or featurecounts
ncores <- args[5]         # MultiCoreParam ncores

if (length(args) == 6) {

    denominator <- args[6]  # denominator for lfc set by user

} else {

    denominator <- ""       # denominator for lfc as default "" meaning is set as first sample condition

}

######################################
######## Define Multicoreparams ######
######################################

BPPARAM <- BiocParallel::MulticoreParam(ncores, stop.on.error = TRUE)

######################################
######## Process sample sheet ########
######################################

# Read in Sample sheet
samps <- read.csv(samplesheet, sep=",", header = TRUE)

# check header of sample sheet
if (!c("sample") %in% colnames(samps) | !c("condition") %in% colnames(samps)) {

    stop("run_dexseq_exon.R Samplesheet must contain 'sample' and 'condition' column headers.", call.=FALSE)

}

# Take only sample and condition columns
samps <- samps[,c("sample", "condition")]

# filter for unique rows based on sample name
samps <- samps[!duplicated(samps[,"sample"]),]

######################################
######## Run DEXseq analysis #########
######################################

# Location of DEXseq count files
countFiles = list.files(countFiles_dir, pattern = ".clean.count.txt$", full.names = TRUE, recursive = TRUE)

# Define models

fullModel <- as.formula("~sample + exon + condition:exon")
reducedModel <- as.formula("~sample + exon")

if (read_method == "htseq"){

    dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(countfiles = countFiles,
                                          sampleData = samps,
                                          design = fullModel,
                                          flattenedfile = flattenedFile)
}

# TODO Implement featurecounts input option
# if (read_method == "featurecounts"){

#     source("load_SubreadOutput.R")

#     dxd <- DEXSeqDataSetFromFeatureCounts(countFiles,
#                                           flattenedfile = flattenedFile,
#                                           sampleData = samps)

# }

dxd <- DEXSeq::estimateSizeFactors(dxd)

dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE, BPPARAM = BPPARAM)

# Looks for condition specific difference in tx proportions
dxd <- DEXSeq::testForDEU(dxd, fullModel = fullModel, reducedModel = reducedModel, BPPARAM = BPPARAM)

# Define sample col used for lfc calculation - at current this is fixed to required "condition" column
fitExpToVar <- "condition"

# Get fold changes based on fitExpToVar col in colData and denominator defines baseline for lfc
dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = fitExpToVar, denominator = denominator, BPPARAM = BPPARAM)

# Get Results
dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)

# Get q vals
qval <- DEXSeq::perGeneQValue(dxr)

# Format q vals
dxr.g <- data.frame(gene = names(qval), qval)

# dxr tsv
keep <- colnames(dxr)[sapply(dxr, class) %in% c("numeric", "character")]
dxr.tsv <- dxr[ , colnames(dxr) %in% keep]

################################
######### Save outputs #########
################################

# dxd
saveRDS(dxd, "dxd_exon.rds")

# results
saveRDS(dxr, "dxr_exon.rds")
write.table(dxr.tsv, "dxr_exon.tsv", sep="\t", quote=FALSE, row.names = TRUE)

# qvals
saveRDS(qval, "qval_exon.rds")
write.table(dxr.g, "dxr_exon.g.tsv", sep="\t", quote=FALSE, row.names = TRUE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DEXSeq")
sessionInfo()
