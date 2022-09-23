#!/usr/bin/env Rscript

library(DEXseq)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 1) {

    stop("Usage: run_dexseq.R <drimseq_filter_rds>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

d = args[1]  # d.rds object from run_drimseq_filter.R

######################################
######## Run DEXseq analysis #########
######################################

# Take pre-filtered sample data from DRIMSeq object
sample.data <- DRIMSeq::samples(d)

# Take count data from same filtered DRIMSeq object
# count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
count.data <- round(counts(d))

# DEXseq made for exon level but modified for DTU
# design specifies "exon" but should be read as "transcript"
# See F1000 workflow for more details:
# https://f1000research.com/articles/7-952

dxd <- DEXSeqDataSet(countData = count.data,
                    sampleData = sample.data,
                    design = ~sample + exon + condition:exon, 
                    featureID = counts(d)$feature_id,
                    groupID = counts(d)$gene_id)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)

# Looks for condition specific difference in tx proportions
dxd <- testForDEU(dxd, reducedModel = ~sample + exon)

# Get Results
dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)

# Format results
columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])

# Get q vals
qval <- perGeneQValue(dxr)

# Format q vals
dxr.g <- data.frame(gene=names(qval),qval)

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
citation("DEXseq")
sessionInfo()