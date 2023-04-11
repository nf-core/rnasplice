#!/usr/bin/env Rscript

library(DEXSeq)

args <- commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 2) {

  stop("Usage: run_dexseq_plot.R <dxr> <n_dexseq_plot>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

dxr <- args[1]
n   <- args[2]

##########################################
######## Plot  and save DEXseq  ##########
##########################################

dxr=readRDS(dxr)
unq <- unique(dxr$groupID)
top_dxr=unique(dxr$groupID[order(dxr$padj)])[1:n]

pdf(file= "dexseq_plot.pdf", width=11.7, height=8.3)

for (i in 1:n)
  plotDEXSeq( dxr, top_dxr[i],legend=TRUE, displayTranscripts=TRUE, splicing=TRUE, cex.axis=1, cex=1, lwd=2)
dev.off()

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DEXSeq")
sessionInfo()
