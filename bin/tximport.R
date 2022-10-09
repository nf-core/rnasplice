#!/usr/bin/env Rscript

library(tximport)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 3) {
  
  stop("Usage: tximport.R <tx2gene> <salmon_out> <sample_name>", call.=FALSE)
  
}

######################################
########### Collect inputs ###########
######################################

tx2gene = args[1]  # "<prefix>_tx2gene.tsv"
path = args[2]
prefix = args[3]

# Read in tx2gene file

if (!file.exists(tx2gene)) {
  
  stop("Usage: tximport.R <tx2gene> <salmon_out> <sample_name> - No tx2gene.tsv specified", call.=FALSE)
  
} else {
  
  # Read in tx2gene file
  rowdata <- read.csv(tx2gene, sep="\t", header = FALSE)
  
  # Set tx2gene header
  colnames(rowdata) <- c("tx", "gene_id", "gene_name")
  
  # Take only first 2 cols
  tx2gene <- rowdata[,1:2]
  
}

# Collect salmon quant files

fns <- list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names <- basename(dirname(fns))
names(fns) <- names

####################################
########### Run Tximport ###########
####################################

# Run Tximport across countsFromAbundance options

txi <- tximport::tximport(fns, type = "salmon", txOut = TRUE,
                          countsFromAbundance = "no")

txi.s <- tximport::tximport(fns, type = "salmon", txOut = TRUE,
                            countsFromAbundance = "scaledTPM")

txi.ls <- tximport::tximport(fns, type = "salmon", txOut = TRUE,
                             countsFromAbundance = "lengthScaledTPM")

####################################################
########### Run Tximport:summarizeToGene ###########
####################################################

# Run summarizeToGene

gi <- tximport::summarizeToGene(txi, tx2gene = tx2gene,
                                countsFromAbundance = "no")

gi.s <- tximport::summarizeToGene(txi, tx2gene = tx2gene,
                                  countsFromAbundance = "scaledTPM")

gi.ls <- tximport::summarizeToGene(txi, tx2gene = tx2gene,
                                   countsFromAbundance="lengthScaledTPM")

####################################################
########### Run Tximport:dtuScaledTPM ##############
####################################################

# Add in tx ids from salmon quants into tx2gene to ensure Tximport:dtuScaledTPM runs

missing_txids <- setdiff(rownames(txi[[1]]),  as.character(tx2gene[["tx"]]))

if (length(missing_txids) > 0) {
  
  message("transcripts missing from tx2gene for Tximport:dtuScaledTPM: ", length(missing_txids))
  
  tx2gene_complete <- rbind(tx2gene, data.frame(tx = missing_txids, gene_id = missing_txids))
  
  tx2gene_complete <- tx2gene_complete[match(rownames(txi[[1]]), as.character(tx2gene_complete[["tx"]])),]
  
  txi.dtu <- tximport::tximport(fns, type = "salmon", tx2gene = tx2gene_complete,
                                txOut = TRUE, countsFromAbundance = "dtuScaledTPM")
} else {
  
  txi.dtu <- tximport::tximport(fns, type = "salmon", tx2gene = tx2gene,
                                txOut = TRUE, countsFromAbundance = "dtuScaledTPM")
}

##############################################################################
####### Check tx2gene tx and txis to ensure consistency prior to output ######
##############################################################################

filter_txi <- function(txi.obj, tx2gene_tsv){
  
  # unpack matrices
  abundanceMatTx <- txi.obj$abundance
  countsMatTx <- txi.obj$counts
  lengthMatTx <- txi.obj$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  if (!any(txId %in% tx2gene_tsv$tx)) {
    txFromFile <- paste0("Example IDs (file): [", paste(head(txId,3),collapse=", "),", ...]")
    txFromTable <- paste0("Example IDs (tx2gene): [", paste(head(tx2gene_tsv$tx,3),collapse=", "),", ...]")
    stop(paste0("
  None of the transcripts in the quantification files are present
  in the first column of tx2gene. Check to see that you are using
  the same annotation for both.\n\n",txFromFile,"\n\n",txFromTable))
  }
  
  # remove transcripts (and genes) not in the rownames of matrices
  tx2gene_tsv <- tx2gene_tsv[tx2gene_tsv$tx %in% txId,]
  ntxmissing <- sum(!txId %in% tx2gene_tsv$tx)
  if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing)
  
  # subset to transcripts in the tx2gene table
  sub.idx <- txId %in% tx2gene_tsv$tx
  abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
  countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
  lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
  
  # resave matrices
  txi.obj$abundance <- abundanceMatTx
  txi.obj$counts <- countsMatTx
  txi.obj$length <- lengthMatTx
  
  return(list(txi.obj, tx2gene_tsv))
}

# Run through filter to ensure txi and tx2gene match for downstream analysis

txi <- filter_txi(txi, tx2gene)[[1]]
txi.s <- filter_txi(txi.s, tx2gene)[[1]]
txi.ls <- filter_txi(txi.ls, tx2gene)[[1]]
txi.dtu <- filter_txi(txi.dtu, tx2gene)[[1]]
tx2gene <- filter_txi(txi, tx2gene)[[2]]

missing_txids <- setdiff(rownames(txi[[1]]),  as.character(tx2gene[["tx"]]))
stopifnot(length(missing_txids) == 0) 

####################################
########### Save Output ############
####################################

# Save out txi rds

saveRDS(txi, paste(c(prefix, "txi.rds"), collapse="."))
saveRDS(txi.s, paste(c(prefix, "txi.s.rds"), collapse="."))
saveRDS(txi.ls, paste(c(prefix, "txi.ls.rds"), collapse="."))
saveRDS(txi.dtu, paste(c(prefix, "txi.dtu.rds"), collapse="."))

# Save out gi rds

saveRDS(gi, paste(c(prefix, "gi.rds"), collapse="."))
saveRDS(gi.s, paste(c(prefix, "gi.s.rds"), collapse="."))
saveRDS(gi.ls, paste(c(prefix, "gi.ls.rds"), collapse="."))

# Save out tpm tsv

# tx level
write.table(cbind.data.frame(tx2gene,txi[["abundance"]]), paste(c(prefix, "transcript_tpm.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,txi[["counts"]]), paste(c(prefix, "transcript_counts.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(tx2gene,txi.s[["abundance"]]), paste(c(prefix, "transcript_tpm_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,txi.s[["counts"]]), paste(c(prefix, "transcript_counts_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(tx2gene,txi.ls[["abundance"]]), paste(c(prefix, "transcript_tpm_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,txi.ls[["counts"]]), paste(c(prefix, "transcript_counts_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(tx2gene,txi.dtu[["abundance"]]), paste(c(prefix, "transcript_tpm_dtu_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,txi.dtu[["counts"]]), paste(c(prefix, "transcript_counts_dtu_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

# gene level
write.table(cbind.data.frame(gene_id = rownames(gi[["abundance"]]), gi[["abundance"]]), paste(c(prefix, "gene_tpm.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(gene_id = rownames(gi[["counts"]]),gi[["counts"]]), paste(c(prefix, "gene_counts.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(gene_id = rownames(gi.s[["abundance"]]),gi.s[["abundance"]]), paste(c(prefix, "gene_tpm_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(gene_id = rownames(gi.s[["counts"]]),gi.s[["counts"]]), paste(c(prefix, "gene_counts_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(gene_id = rownames(gi.ls[["abundance"]]),gi.ls[["abundance"]]), paste(c(prefix, "gene_tpm_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(gene_id = rownames(gi.ls[["counts"]]),gi.ls[["counts"]]), paste(c(prefix, "gene_counts_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

# output new tx2gene
write.table(tx2gene, "tximport.tx2gene.tsv", sep="\t", quote=FALSE, row.names = FALSE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("tximport")
sessionInfo()