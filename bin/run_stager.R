#!/usr/bin/env Rscript

library(stageR)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 2) {

    stop("Usage: run_stager.R <dexseq_results_rds/drimseq_results_rds> <analysis_type> <qvals>|<res.txp>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

results = args[1]
analysis_type = args[2]  

#########################################
####### Run stageR postprocessing #######
#########################################

strp <- function(x) substr(x,1,15)

if (analysis_type == "dexseq"){

    dxr <- readRDS(results)
    
    qval <- args[3]

    qval <- readRDS(qval)

    dxr <- as.data.frame(dxr)

    pConfirmation <- matrix(dxr$pvalue,ncol=1)
    dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")

    pScreen <- qval
    names(pScreen) <- strp(names(pScreen))

    tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
    for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

    res_pval <- as.data.frame(dxr[,1:7])

} else if (analysis_type == "drimseq"){

    res <- readRDS(results)

    res.txp <- args[3]

    res.txp <- readRDS(res.txp)

    res <- as.data.frame(res)

    pConfirmation <- matrix(res.txp$pvalue, ncol=1)
    rownames(pConfirmation) <- strp(res.txp$feature_id)

    pScreen <- res$pvalue
    names(pScreen) <- strp(res$gene_id)

    tx2gene <- res.txp[,c("feature_id", "gene_id")]
    for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

    res_pval <- as.data.frame(res[,1:6])
}

# Run StageR 

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)

stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=TRUE)

suppressWarnings({
    stageR.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                      onlySignificantGenes=FALSE)
})

################################
######### Save outputs #########
################################

# StageR object
saveRDS(stageRObj, paste0(analysis_type,".stageRObj.rds"))

# Adjusted P values
saveRDS(stageR.padj, paste0(analysis_type,".stageR.padj.rds"))
write.table(stageR.padj, paste0(analysis_type,".stageR.padj.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

# Combine outputs of stager and drim/dexseq
colnames(stageR.padj) <- paste("stageR", colnames(stageR.padj), sep = ".")
colnames(res_pval) <- paste(analysis_type, colnames(res_pval), sep = ".")

write.table(cbind.data.frame(res_pval,stageR.padj), paste0(analysis_type,".combined.stageR.padj.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("stageR")
sessionInfo()