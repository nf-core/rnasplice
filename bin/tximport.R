#!/usr/bin/env Rscript

library(tximport)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 2) {
  
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
####### Match tx2gene txids and salmon quants ######
####################################################

# Ensure all tx ids from salmon quants present in tx2gene

missing_txids <- setdiff(rownames(txi[[1]]),  as.character(tx2gene[["tx"]]))

if (length(missing_txids) > 0) {
  
  tx2gene <- rbind(tx2gene, data.frame(tx = missing_txids, gene_id = missing_txids))
}

tx2gene <- tx2gene[match(rownames(txi[[1]]), as.character(tx2gene[["tx"]])),]

####################################################
########### Run Tximport:dtuScaledTPM ##############
####################################################

# Run Tximport across countsFromAbundance options

txi.dtu <- tximport::tximport(fns, type = "salmon", tx2gene = tx2gene,
                              txOut = TRUE, countsFromAbundance = "dtuScaledTPM")

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
write.table(cbind.data.frame(tx2gene,gi[["abundance"]]), paste(c(prefix, "gene_tpm.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,gi[["counts"]]), paste(c(prefix, "gene_counts.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(tx2gene,gi.s[["abundance"]]), paste(c(prefix, "gene_tpm_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,gi.s[["counts"]]), paste(c(prefix, "gene_counts_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

write.table(cbind.data.frame(tx2gene,gi.ls[["abundance"]]), paste(c(prefix, "gene_tpm_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(cbind.data.frame(tx2gene,gi.ls[["counts"]]), paste(c(prefix, "gene_counts_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

# output new tx2gene
write.table(tx2gene, "tximport.tx2gene.tsv", sep="\t", quote=FALSE, row.names = FALSE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("tximport")
sessionInfo()