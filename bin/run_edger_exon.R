#!/usr/bin/env Rscript


# Read command arguments

argv <- commandArgs(trailingOnly = TRUE)

argc <- length(argv)


# Parse command arguments

featurecounts <- argv[1]

samplesheet <- argv[2]

contrastsheet <- argv[3]

ntop  <- argv[4]


# Load required packages

library(edgeR)


# Read samples table

samples <- read.csv(samplesheet)

samples <- samples[, c("sample", "condition"), drop = FALSE]

samples <- unique(samples)


# Read contrasts table

contrasts <- read.csv(contrastsheet)

contrasts <- contrasts[, c("contrast", "treatment", "control"), drop = FALSE]


# Read featureCounts files

files <- paste0(featurecounts, "/", samples$sample, ".featureCounts.txt")

data <- lapply(files, read.delim, comment.char = "#")

data <- Reduce(merge, data)


# Extract counts matrix

counts <- data[, -(1:6)]

counts <- as.matrix(counts)

colnames(counts) <- samples$sample


# Extract genes annotation

genes <- data[, 1:6]


# Create DGEList object

DGEList <- DGEList(
    counts  = counts,
    samples = samples,
    group   = samples$condition,
    genes   = genes,
)


# Normalization

keep <- filterByExpr(DGEList, group = DGEList$samples$group)

DGEList <- DGEList[keep, , keep.lib.sizes = FALSE]

DGEList <- calcNormFactors(DGEList)


# Create design matrix

group <- factor(DGEList$samples$group)

design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)


# Create contrasts matrix

groups <- levels(group)

names <- contrasts$contrast

contrasts <- data.frame(A = contrasts$treatment, B = contrasts$control)

contrasts <- apply(contrasts, 1, paste, collapse = "-")

contrasts <- makeContrasts(contrasts = contrasts, levels = groups)

colnames(contrasts) <- names

contrasts <- contrasts[, colSums(contrasts != 0) > 0]


# Estimate dispersions by empirical Bayes

DGEList <- estimateDisp(DGEList, design)


# Fit log-linear model to count data

DGEGLM <- glmQLFit(DGEList, design)


# Test for differential exon expression

DGELRT.exprs <- mapply(
    FUN      = glmQLFTest,
    contrast = asplit(contrasts, MARGIN = 2),
    MoreArgs = list(glmfit = DGEGLM),
    SIMPLIFY = FALSE
)

results.exprs <- lapply(DGELRT.exprs, topTags, n = Inf, sort.by = "none")


# Test for differential exon usage

DGELRT.usage <- mapply(
    FUN      = diffSpliceDGE,
    contrast = asplit(contrasts, MARGIN = 2),
    MoreArgs = list(glmfit = DGEGLM, geneid = "Geneid", exonid = "Start"),
    SIMPLIFY = FALSE
)

results.usage <- list(
    simes = lapply(DGELRT.usage, topSpliceDGE, test = "Simes", number = Inf),
    gene  = lapply(DGELRT.usage, topSpliceDGE, test = "gene", number = Inf),
    exon  = lapply(DGELRT.usage, topSpliceDGE, test = "exon", number = Inf)
)


# Save objects to disk

saveRDS(DGEList, file = "DGEList.rds")

saveRDS(DGEGLM,  file = "DGEGLM.rds")

saveRDS(DGELRT.exprs,  file = "DGELRT.exprs.rds")

saveRDS(DGELRT.usage,  file = "DGELRT.usage.rds")


# Save results to disk

mapply(
    write.csv,
    x = results.exprs,
    file = paste0("contrast_", colnames(contrasts), ".exprs.csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
    write.csv,
    x = results.usage$simes,
    file = paste0("contrast_", colnames(contrasts), ".usage.simes.csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
    write.csv,
    x = results.usage$gene,
    file = paste0("contrast_", colnames(contrasts), ".usage.gene.csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
    write.csv,
    x = results.usage$exon,
    file = paste0("contrast_", colnames(contrasts), ".usage.exon.csv"),
    MoreArgs = list(quote = FALSE, row.names = FALSE)
)


# Save plots to disk

write.plotSpliceDGE <- function(results, file, lrt, n = 10) {

    ids <- head(results$Geneid, n = n)

    pdf(file = file, width = 11.7, height = 8.3)

    for (i in ids) { plotSpliceDGE(lrt, geneid = i) }

    dev.off()

}

mapply(
    write.plotSpliceDGE,
    results = results.usage$simes,
    file = paste0("contrast_", colnames(contrasts), ".usage.simes.pdf"),
    lrt = DGELRT.usage,
    MoreArgs = list(n = ntop)
)

mapply(
    write.plotSpliceDGE,
    results = results.usage$gene,
    file = paste0("contrast_", colnames(contrasts), ".usage.gene.pdf"),
    lrt = DGELRT.usage,
    MoreArgs = list(n = ntop)
)

mapply(
    write.plotSpliceDGE,
    results = results.usage$exon,
    file = paste0("contrast_", colnames(contrasts), ".usage.exon.pdf"),
    lrt = DGELRT.usage,
    MoreArgs = list(n = ntop)
)


# Print session information

citation("edgeR")

sessionInfo()
