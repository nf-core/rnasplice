#!/usr/bin/env Rscript

library(stageR)

args = commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 2) {

    stop("Usage: run_stager.R <dexseq_results_rds/drimseq_results_rds> <analysis_type> <tx2gene> <qvals>|<res.txp>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

results = args[1]
analysis_type = args[2]  
tx2gene = args[3]

#########################################
####### Run stageR postprocessing #######
#########################################

strp <- function(x) substr(x,1,15)

if (analysis_type == "dexseq"){

    dxr <- readRDS(results)
    
    qval <- args[4]

    qval <- readRDS(qval)

    pConfirmation <- matrix(dxr$pvalue,ncol=1)
    dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")

    pScreen <- qval
    names(pScreen) <- strp(names(pScreen))

    tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
    for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

} else if (analysis_type == "drimseq"){

    res <- readRDS(results)

    res.txp <- args[4]
    
    res.txp <- readRDS(res.txp)

    pConfirmation <- matrix(res.txp$pvalue, ncol=1)
    rownames(pConfirmation) <- strp(res.txp$feature_id)

    pScreen <- res$pvalue
    names(pScreen) <- strp(res$gene_id)

    tx2gene <- res.txp[,c("feature_id", "gene_id")]
    for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
}

# Run StageR 

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=FALSE, tx2gene=tx2gene)

stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)

suppressWarnings({
    stageR.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                    onlySignificantGenes=TRUE)
})

################################
######### Save outputs #########
################################

# StageR object
saveRDS(stageRObj, paste0(analysis_type,".stageRObj.rds"))

# Adjusted P values
saveRDS(stageR.padj, paste0(analysis_type,".stageR.padj.rds"))
write.table(stageR.padj, paste0(analysis_type,".stageR.padj.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("stageR")
sessionInfo()