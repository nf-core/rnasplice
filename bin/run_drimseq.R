#!/usr/bin/env Rscript

library(DRIMSeq)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 3) {

    stop("Usage: run_drimseq.R <drimseq_filter_rds> <samplesheet> <coef>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

d = args[1]  # d.rds object from run_drimseq_filter.R
samplesheet = args[2] # samplesheet
coef = args[3] # coefficient for comparison

######################################
########### Run DRIMSeq analysis #####
######################################

## Specify model
design_full <- model.matrix(~condition, data = DRIMSeq::samples(d))

# Run DRIMseq
d <- dmPrecision(d, design = design_full)
d <- dmFit(d, design = design_full)
d <- dmTest(d, coef = coef)

# Get results
res <- DRIMSeq::results(d)

# feature level
res.txp <- DRIMSeq::results(d, level="feature")

# Filter (remove NAs)
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)


# idx <- which(res$adj_pvalue < 0.05)[1]
# p <- plotProportions(d, res$gene_id[idx], "condition")

# pdf("Example_DTU.pdf")
# plot(p)
# dev.off()

################################
######### Save outputs #########
################################

saveRDS(d, "drimseq_d.rds")
saveRDS(res, "res.rds")
saveRDS(res.txp, "res.txp.rds")

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DRIMSeq")
sessionInfo()
