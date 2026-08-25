#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

library(tximport)

######################################
## Functions                        ##
######################################

parse_args <- function(x) {
  args_list <- unlist(strsplit(x, " ?--")[[1]])[-1]
  args_vals <- lapply(
    args_list, function(y) scan(text = y, what = "character", quiet = TRUE)
)

  # Ensure the option vectors are length 2 (key/value) to catch empty ones
  args_vals <- lapply(args_vals, function(z) {
    length(z) <- 2
    z
  })

  parsed_args <- structure(
    lapply(
      args_vals, function(y) y[2]), names = lapply(args_vals, function(y) y[1]
    )
  )
  parsed_args[!is.na(parsed_args)]
}

######################################
########### Collect inputs ###########
######################################

# Variables below are interpolated by Nextflow when rendering this template.

tx2gene <- "$tx2gene" # "<prefix>_tx2gene.tsv"
path <- "salmon"
prefix <- ""
if ("$task.ext.prefix" != "null") {
  prefix <- "$task.ext.prefix"
} else if ("$meta.id" != "null") {
  prefix <- "$meta.id"
}

args_opt <- parse_args("$task.ext.args")
if (!is.null(args_opt[["ignore_tx_version"]])) {
  ignore_tx_version <- as.logical(args_opt[["ignore_tx_version"]])
} else {
  ignore_tx_version <- NA
}

if (is.na(ignore_tx_version)) {
  message(
    "No --ignore_tx_version argument specified. Defaulting to FALSE."
  )
  ignore_tx_version <- FALSE
}

# Read in tx2gene file

if (!file.exists(tx2gene)) {
  stop("No tx2gene.tsv specified", call. = FALSE)
} else {
  # Read in tx2gene file
  rowdata <- read.csv(tx2gene, sep = "\t", header = FALSE)

  # Set tx2gene header
  colnames(rowdata) <- c("tx", "gene_id", "gene_name")

  # Take only first 2 cols
  tx2gene <- rowdata[, 1:2]
}

# If ignore_tx_version is TRUE, remove version numbers from tx ids in tx2gene
# and salmon quant files
if (ignore_tx_version) {
  message(
    "Ignoring transcript version numbers in tx2gene and salmon quant files."
  )

  # Remove version numbers from tx2gene
  tx2gene\$tx <- sub("[.][0-9]+\$", "", tx2gene\$tx)
}

# Collect salmon quant files

fns <- list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names <- basename(dirname(fns))
names(fns) <- names

####################################
########### Run Tximport ###########
####################################

# Run Tximport across countsFromAbundance options

txi <- tximport::tximport(fns,
  type = "salmon", txOut = TRUE,
  countsFromAbundance = "no", ignoreTxVersion = ignore_tx_version
)

txi.s <- tximport::tximport(fns,
  type = "salmon", txOut = TRUE,
  countsFromAbundance = "scaledTPM", ignoreTxVersion = ignore_tx_version
)

txi.ls <- tximport::tximport(fns,
  type = "salmon", txOut = TRUE,
  countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = ignore_tx_version
)

####################################################
########### Run Tximport:summarizeToGene ###########
####################################################

# Run summarizeToGene

gi <- tximport::summarizeToGene(txi,
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = ignore_tx_version
)

gi.s <- tximport::summarizeToGene(txi,
  tx2gene = tx2gene,
  countsFromAbundance = "scaledTPM",
  ignoreTxVersion = ignore_tx_version
)

gi.ls <- tximport::summarizeToGene(txi,
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = ignore_tx_version
)

####################################################
########### Run Tximport:dtuScaledTPM ##############
####################################################

# Add in tx ids from salmon quants into tx2gene to ensure Tximport:dtuScaledTPM
# runs

txi_rownames <- if (ignore_tx_version) {
  sub("[.][0-9]+\$", "", rownames(txi[[1]]))
} else {
  rownames(txi[[1]])
}

missing_txids <- setdiff(txi_rownames, as.character(tx2gene[["tx"]]))

if (length(missing_txids) > 0) {
  message(
    "transcripts missing from tx2gene for Tximport:dtuScaledTPM: ",
    length(missing_txids)
  )

  tx2gene_complete <- rbind(
    tx2gene, data.frame(tx = missing_txids, gene_id = missing_txids)
  )

  tx2gene_complete <- tx2gene_complete[
    match(rownames(txi[[1]]), as.character(tx2gene_complete[["tx"]])),
  ]

  txi.dtu <- tximport::tximport(
    fns,
    type = "salmon", tx2gene = tx2gene_complete,
    txOut = TRUE, countsFromAbundance = "dtuScaledTPM",
    ignoreTxVersion = ignore_tx_version
  )

} else {
  txi.dtu <- tximport::tximport(fns,
    type = "salmon", tx2gene = tx2gene,
    txOut = TRUE, countsFromAbundance = "dtuScaledTPM",
    ignoreTxVersion = ignore_tx_version
  )
}

##############################################################################
####### Check tx2gene tx and txis to ensure consistency prior to output ######
##############################################################################

### The below function has been taken and modified from the tximport package
### and summarizeToGene function
### citation("tximport")

filter_txi <- function(txi.obj, tx2gene_tsv) {
  # unpack matrices
  abundanceMatTx <- txi.obj\$abundance
  countsMatTx <- txi.obj\$counts
  lengthMatTx <- txi.obj\$length

  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))

  if (!any(txId %in% tx2gene_tsv\$tx)) {
    txFromFile <- paste0(
      "Example IDs (file): [", paste(head(txId, 3), collapse = ", "), ", ...]"
    )
    txFromTable <- paste0("Example IDs (tx2gene): [",
      paste(head(tx2gene_tsv\$tx,3),collapse=", "),", ...]"
    )

    stop(
      paste0("
        None of the transcripts in the quantification files are present
        in the first column of tx2gene. Check to see that you are using
        the same annotation for both.\n\n", txFromFile, "\n\n", txFromTable
      )
    )
  }

  # remove transcripts (and genes) not in the rownames of matrices
  ntx2genemissing <- sum(!tx2gene_tsv\$tx %in% txId)
  if (ntx2genemissing > 0) {
    message("Filtering transcripts from tx2gene: ", ntx2genemissing)
  }

  tx2gene_tsv <- tx2gene_tsv[tx2gene_tsv\$tx %in% txId,]

  # subset to transcripts in the tx2gene table
  ntxmissing <- sum(!txId %in% tx2gene_tsv\$tx)
  if (ntxmissing > 0) {
    message("Filtering transcripts from txi: ", ntxmissing)
  }

  sub.idx <- txId %in% tx2gene_tsv\$tx
  abundanceMatTx <- abundanceMatTx[sub.idx, , drop = FALSE]
  countsMatTx <- countsMatTx[sub.idx, , drop = FALSE]
  lengthMatTx <- lengthMatTx[sub.idx, , drop = FALSE]

  # resave matrices
  txi.obj\$abundance <- abundanceMatTx
  txi.obj\$counts <- countsMatTx
  txi.obj\$length <- lengthMatTx

  tx2gene_tsv <- tx2gene_tsv[match(rownames(txi.obj\$abundance),
    as.character(tx2gene_tsv[["tx"]])),
  ]

  if (!all(rownames(txi.obj\$abundance) == tx2gene_tsv\$tx)) {
    stop("Stop: tx2gene and txi rownames do not match - Cannot proceed with merge.")
  }

  return(list(txi.obj, tx2gene_tsv))
}

# Run through filter to ensure txi and tx2gene match for downstream analysis

txi <- filter_txi(txi, tx2gene)[[1]]
txi.s <- filter_txi(txi.s, tx2gene)[[1]]
txi.ls <- filter_txi(txi.ls, tx2gene)[[1]]
txi.dtu <- filter_txi(txi.dtu, tx2gene)[[1]]
tx2gene <- filter_txi(txi, tx2gene)[[2]]

missing_txids <- setdiff(rownames(txi[[1]]), as.character(tx2gene[["tx"]]))
stopifnot(length(missing_txids) == 0)

####################################
########### Save Output ############
####################################

# Save out txi rds

saveRDS(txi, paste(c(prefix, "txi.rds"), collapse = "."))
saveRDS(txi.s, paste(c(prefix, "txi.s.rds"), collapse = "."))
saveRDS(txi.ls, paste(c(prefix, "txi.ls.rds"), collapse = "."))
saveRDS(txi.dtu, paste(c(prefix, "txi.dtu.rds"), collapse = "."))

# Save out gi rds

saveRDS(gi, paste(c(prefix, "gi.rds"), collapse = "."))
saveRDS(gi.s, paste(c(prefix, "gi.s.rds"), collapse = "."))
saveRDS(gi.ls, paste(c(prefix, "gi.ls.rds"), collapse = "."))

# Save out tpm tsv

# tx level #####################################################################

write.table(
  cbind.data.frame(tx2gene, txi[["abundance"]]),
  paste(c(prefix, "transcript_tpm.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cbind.data.frame(tx2gene, txi[["counts"]]),
  paste(c(prefix, "transcript_counts.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cbind.data.frame(tx2gene, txi.s[["abundance"]]),
  paste(c(prefix, "transcript_tpm_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(cbind.data.frame(tx2gene, txi.s[["counts"]]),
  paste(c(prefix, "transcript_counts_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(cbind.data.frame(tx2gene, txi.ls[["abundance"]]),
  paste(c(prefix, "transcript_tpm_length_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(cbind.data.frame(tx2gene, txi.ls[["counts"]]),
  paste(c(prefix, "transcript_counts_length_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(cbind.data.frame(tx2gene, txi.dtu[["abundance"]]),
  paste(c(prefix, "transcript_tpm_dtu_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(cbind.data.frame(tx2gene, txi.dtu[["counts"]]),
  paste(c(prefix, "transcript_counts_dtu_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# gene level ###################################################################

write.table(
  cbind.data.frame(gene_id = rownames(gi[["abundance"]]), gi[["abundance"]]),
  paste(c(prefix, "gene_tpm.tsv"), collapse = "."),
  sep = "\t", quote = FALSE,
  row.names = FALSE
)

write.table(
  cbind.data.frame(gene_id = rownames(gi[["counts"]]), gi[["counts"]]),
  paste(c(prefix, "gene_counts.tsv"), collapse = "."),
  sep = "\t",
  quote = FALSE, row.names = FALSE
)

write.table(
  cbind.data.frame(
    gene_id = rownames(gi.s[["abundance"]]),
    gi.s[["abundance"]]
  ),
  paste(c(prefix, "gene_tpm_scaled.tsv"), collapse = "."),
  sep = "\t",
  quote = FALSE, row.names = FALSE
)

write.table(
  cbind.data.frame(
    gene_id = rownames(gi.s[["counts"]]), gi.s[["counts"]]
  ),
  paste(c(prefix, "gene_counts_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cbind.data.frame(
    gene_id = rownames(gi.ls[["abundance"]]),
    gi.ls[["abundance"]]
  ), paste(c(prefix, "gene_tpm_length_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cbind.data.frame(gene_id = rownames(gi.ls[["counts"]]), gi.ls[["counts"]]),
  paste(c(prefix, "gene_counts_length_scaled.tsv"), collapse = "."),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# output new tx2gene
write.table(
  tx2gene, "tximport.tx2gene.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# output single tpm tsv for suppa downstream
suppa_tpm <- txi[["abundance"]]
write.table(
  suppa_tpm, "suppa_tpm.txt",
  sep = "\t", quote = FALSE, row.names = TRUE
)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("tximport")
sessionInfo()

####################################
########### Versions ###############
####################################

writeLines(
  c(
    '"${task.process}":',
    # paste("    r-base:", paste0(R.version\$major, ".", R.version\$minor)),
    paste("    bioconductor-tximeta:", as.character(packageVersion("tximeta"))),
    paste(
      "    bioconductor-tximport:",
      as.character(packageVersion("tximport"))
    )
  ),
  "versions.yml"
)
