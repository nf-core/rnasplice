#!/usr/bin/env Rscript

# Parse command arguments

argv <- commandArgs(trailingOnly = TRUE)

argc <- length(argv)

featurecounts <- argv[1]

samplesheet <- argv[2]

# Load required packages

library(edgeR)

# Read samples table

samples <- read.csv(samplesheet)

samples <- samples[, c("sample", "condition"), drop = FALSE]

samples <- unique(samples)

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

DGEList <- calcNormFactors(DGEList)

# Create design matrix

group <- factor(DGEList$samples$group)

design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)

# Create contrast matrix

groups <- levels(group)

combinations <- expand.grid(A = groups, B = groups)

contrasts <- apply(combinations, 1, paste, collapse = "-")

contrasts <- makeContrasts(contrasts = contrasts, levels = groups)

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
  file = paste0("contrast_", names(results.exprs), ".exprs.csv"),
  MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
  write.csv,
  x = results.usage$simes,
  file = paste0("contrast_", names(results.usage$simes), ".usage.simes.csv"),
  MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
  write.csv,
  x = results.usage$gene,
  file = paste0("contrast_", names(results.usage$gene), ".usage.gene.csv"),
  MoreArgs = list(quote = FALSE, row.names = FALSE)
)

mapply(
  write.csv,
  x = results.usage$exon,
  file = paste0("contrast_", names(results.usage$exon), ".usage.exon.csv"),
  MoreArgs = list(quote = FALSE, row.names = FALSE)
)
