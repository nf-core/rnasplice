#!/usr/bin/env Rscript

library(DRIMSeq)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 3) {
  
  stop("Usage: run_drimseq_filter.R <txi> <tx2gene> <samplesheet>", call.=FALSE)
  
}

######################################
########### Collect inputs ###########
######################################

txi = args[1]  # txi.rds object from tximport
tx2gene = args[2] # tx2gene information
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
rowdata <- read.csv(tx2gene, sep="\t", header = FALSE)

# Set tx2gene header
colnames(rowdata) <- c("tx", "gene_id", "gene_name")

# Take only first 2 cols - tx and gene_id
tx2gene <- rowdata[,1:2]

######################################
######## Process sample sheet ########
######################################

# Read in Sample sheet
samps <- read.csv(samplesheet, sep=",", header = TRUE)

# check header of sample sheet
if (!c("sample", "condition") %in% colnames(samps)) {
  
  stop("run_drimseq_filter.R Samplesheet must contain 'sample' and 'condition' column headers.", call.=FALSE)
  
}

# Take only sample and condition columns
samps <- samps[,c("sample", "condition")]

# filter for unique rows based on sample name
samps <- samps[,!duplicated(samps[,"sample"])]

######################################
#### Get Counts from tximport txi ####
######################################

# Load in txi from tximport module
txi <- readRDS(txi)

# Take the counts from txi (will be scaledTPM or dtuScaledTPM)
cts <- txi$counts

# Filter for txs with > 0 counts across all samples
cts <- cts[rowSums(cts) > 0,]

# Create counts data frame used downstream
counts <- data.frame(gene_id = tx2gene$gene_id,
                      feature_id = tx2gene$tx,
                      cts)

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