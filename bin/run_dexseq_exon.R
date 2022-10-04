#!/usr/bin/env Rscript

library(DEXSeq)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 1) {
  
  stop("Usage: run_dexseq_exon.R <countFiles> <flattenedFile> <samplesheet> <read_method>", call.=FALSE)
  
}

######################################
########### Collect inputs ###########
######################################

countFiles_dir <- args[1] # count files
flattenedFile <- args[2]  # gff
samplesheet <- args[3]    # samplesheet
read_method <- args[4]    # either HTSeq or featurecounts

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
fitExpToVar <- "condition"
denominator <- ""

if (read_method == "htseq"){
    
    dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(countfiles = countFiles,
                                          sampleData = samps,
                                          design = fullModel,
                                          flattenedfile = flattenedFile)
}

dxd <- DEXSeq::estimateSizeFactors(dxd)

dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE)

# Looks for condition specific difference in tx proportions
dxd <- DEXSeq::testForDEU(dxd, fullModel = fullModel, reducedModel = reducedModel)

# Get fold changes based on condition col in colData
dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = fitExpToVar, denominator = denominator)

# Get Results
dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)

# Get q vals
qval <- DEXSeq::perGeneQValue(dxr)

# Format q vals
dxr.g <- data.frame(gene = names(qval), qval)

# dxr tsv
drop <- c("genomicData","countData")
dxr.tsv <- as.data.frame(dxr)[ , !(names(as.data.frame(dxr)) %in% drop)]

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