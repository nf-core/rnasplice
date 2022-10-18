txt <- c(
  "RAP1_IAA_30M_REP1.featureCounts.txt",
  "RAP1_UNINDUCED_REP1.featureCounts.txt",
  "RAP1_UNINDUCED_REP2.featureCounts.txt",
  "WT_REP1.featureCounts.txt",
  "WT_REP2.featureCounts.txt"
)

dat <- lapply(txt, read.delim, comment.char = "#")

out <- Reduce(merge, dat)

GENES <- seq(1, 6)
COUNTS <- seq(7, ncol(out))

genes <- out[, GENES]
counts <- as.matrix(out[, COUNTS])

edgeR(
  counts = as.matrix(counts),
  genes = genes
)