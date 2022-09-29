#!/usr/bin/env Rscript

library(DRIMSeq)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 3) {
  
  stop("Usage: run_drimseq_filter.R <txi> <tximport_tx2gene> <samplesheet>", call.=FALSE)
  
}

######################################
########### Collect inputs ###########
######################################

txi = args[1]  # txi.rds object from tximport
tximport_tx2gene = args[2] # tx2gene information
samplesheet = args[3] # samplesheet

# Filter params
min_samps_gene_expr = args[4] # Default: 12
min_samps_feature_expr = args[5] # Default: 6
min_samps_feature_prop = args[6] # Default: 6
min_feature_expr = args[7] # Default: 10
min_feature_prop = args[8] # Default: 0.1
min_gene_expr = args[9] # Default: 10

######################################
########## Process tx2gene ###########
######################################

# Read in tx2gene file
tx2gene <- read.csv(tximport_tx2gene, sep="\t", header = TRUE)

# # Set tx2gene header
# colnames(rowdata) <- c("tx", "gene_id", "gene_name")

# # Take only first 2 cols - tx and gene_id
# tximport_tx2gene <- rowdata[,1:2]

######################################
######## Process sample sheet ########
######################################

# Read in Sample sheet
samps <- read.csv(samplesheet, sep=",", header = TRUE)

# check header of sample sheet
if (!c("sample") %in% colnames(samps) | !c("condition") %in% colnames(samps)) {
  
  stop("run_drimseq_filter.R Samplesheet must contain 'sample' and 'condition' column headers.", call.=FALSE)
  
}

# Take only sample and condition columns
samps <- samps[,c("sample", "condition")]

# filter for unique rows based on sample name
samps <- samps[!duplicated(samps[,"sample"]),]

# Change name of cols for DRIMseq 
colnames(samps) <- c("sample_id", "condition")

######################################
#### Get Counts from tximport txi ####
######################################

# Load in txi from tximport module
txi <- readRDS(txi)

# Take the counts from txi (will be scaledTPM or dtuScaledTPM)
cts <- txi$counts

# ensure tx2gene and txi cts match
tx2gene <- tx2gene[match(rownames(cts),tx2gene$tx),]

# check header of sample sheet
if (!all(rownames(cts) == tx2gene$tx)) {
  
  stop("run_drimseq_filter.R Tx2gene rownames and txi rownames must match.", call.=FALSE)
  
}

# Create counts data frame used downstream
counts <- data.frame(gene_id = tx2gene$gene_id,
                      feature_id = tx2gene$tx,
                      cts)


# Filter for txs with > 0 counts across all samples
counts <- counts[rowSums(counts[,(3:ncol(counts))]) > 0,]

######################################
########### Run DRIMSeq filter #######
######################################

d <- DRIMSeq::dmDSdata(counts = counts, samples = samps)

d <- DRIMSeq::dmFilter(d,
                min_samps_feature_expr = min_samps_feature_expr, 
                min_feature_expr = min_feature_expr,
                min_samps_feature_prop = min_samps_feature_prop, 
                min_feature_prop = min_feature_prop,
                min_samps_gene_expr = min_samps_gene_expr, 
                min_gene_expr = min_gene_expr)

################################
######### Save outputs #########
################################

saveRDS(d, "d.rds")

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DRIMSeq")
sessionInfo()