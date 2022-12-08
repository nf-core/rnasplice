#!/usr/bin/env Rscript

library(DEXSeq)
library(BiocParallel)

args <- commandArgs(trailingOnly=TRUE)

# Check args provided

if (length(args) < 6) {

    stop("Usage: run_dexseq_exon.R <countFiles_dir> <gff> <gtf> <samplesheet> <read_method> <ncores> <denominator>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

countFiles_dir     <- args[1] # count files
flattenedFile_gff  <- args[2] # gff
flattenedFile_gtf  <- args[3] # gtf
samplesheet        <- args[4] # samplesheet
read_method        <- args[5] # either HTSeq or featurecounts
ncores             <- args[6] # MultiCoreParam ncores

if (length(args) == 7) {

    denominator <- args[7]  # denominator for lfc set by user

} else {

    denominator <- ""       # denominator for lfc as default "" meaning is set as first sample condition

}

######################################
######## Define Multicoreparams ######
######################################

BPPARAM <- BiocParallel::MulticoreParam(ncores, stop.on.error = TRUE)

######################################
######## Process sample sheet ########
######################################

# Read in Sample sheet
samps <- read.csv(samplesheet, sep=",", header = TRUE)

# check header of sample sheet
if (!c("sample") %in% colnames(samps) | !c("condition") %in% colnames(samps)) {

    stop("run_dexseq_exon.R Samplesheet must contain 'sample' and 'condition' column headers.", call.=FALSE)

}

# Take only sample and condition columns
samps <- samps[,c("sample", "condition")]

# filter for unique rows based on sample name
samps <- samps[!duplicated(samps[,"sample"]),]

#################################################
######## Define featurescounts function #########
#################################################

## Load Fcount output from : DEXSeq_after_Fcount.R into DEXSeq
## Copyright 2015 Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de). Licence: GPLv3.

## Read Fcount output and convert to dxd
DEXSeqDataSetFromFeatureCounts <- function (countfiles, sampleData,
                                            design = ~sample + exon + condition:exon, flattenedfile = NULL)

    {
        # Take a fcount file and convert it to dcounts for dexseq
        message("Reading and adding Exon IDs for DEXSeq")
        for (file in countfiles){
            if (!exists("dcounts")){
                dcounts <- read.table(file, skip=2)
                dcounts <- dcounts[,c(1,7)]
                } else {
                    temp_dcounts <-read.table(file, skip=2)
                    temp_dcounts <- temp_dcounts[,c(1,7)]
                    dcounts <- cbind(dcounts, temp_dcounts[2])
                    rm(temp_dcounts)
                }
            }
        id <- as.character(dcounts[,1])
        n <- id
        split(n,id) <- lapply(split(n ,id), seq_along )
        rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
        dcounts <- dcounts[,2:ncol(dcounts)]

        #dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name

        ## get genes and exon names out
        splitted <- strsplit(rownames(dcounts), ":")
        exons <- sapply(splitted, "[[", 2)
        genesrle <- sapply(splitted, "[[", 1)

        ## parse the flattened file
        if (!is.null(flattenedfile)) {
            aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
            colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
            aggregates$strand <- gsub("\\.", "*", aggregates$strand)
            aggregates <- aggregates[which(aggregates$class == "exon"), ]
            aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
            aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
            # trim the gene_ids to 255 chars in order to match with featurecounts
            longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
            warning(paste0(longIDs, " aggregate geneIDs were found truncated in featureCounts output"), call. = FALSE)
            aggregates$gene_id <- substr(aggregates$gene_id,1,255)

            transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
            transcripts <- strsplit(transcripts, "\\+")
            exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", aggregates$attr)
            exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
            names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")

            names(transcripts) <- names(exoninfo)
            if (!all(rownames(dcounts) %in% names(exoninfo))) {
                    stop("Count files do not correspond to the flattened annotation file")
            }
            matching <- match(rownames(dcounts), names(exoninfo))
            stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
            stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
            dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, genesrle, exoninfo[matching], transcripts[matching])
            return(dxd)
        }
        else {
            dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, genesrle)
            return(dxd)
        }

}

######################################
######## Run DEXseq analysis #########
######################################

# Location of DEXseq count files
countFiles = list.files(countFiles_dir, pattern = ".txt$", full.names = TRUE, recursive = TRUE)

# Define models

fullModel <- as.formula("~sample + exon + condition:exon")
reducedModel <- as.formula("~sample + exon")

if (read_method == "htseq"){

    dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(countfiles = countFiles,
                                        sampleData = samps,
                                        design = fullModel,
                                        flattenedfile = flattenedFile_gff)
}

if (read_method == "featurecounts"){

    dxd <- DEXSeqDataSetFromFeatureCounts(countfile = countFiles,
                                        sampleData = samps,
                                        design = fullModel,
                                        flattenedfile = flattenedFile_gtf)
}

dxd <- DEXSeq::estimateSizeFactors(dxd)

dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE, BPPARAM = BPPARAM)

# Looks for condition specific difference in tx proportions
dxd <- DEXSeq::testForDEU(dxd, fullModel = fullModel, reducedModel = reducedModel, BPPARAM = BPPARAM)

# Define sample col used for lfc calculation - at current this is fixed to required "condition" column
fitExpToVar <- "condition"

# Get fold changes based on fitExpToVar col in colData and denominator defines baseline for lfc
dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = fitExpToVar, denominator = denominator, BPPARAM = BPPARAM)

# Get Results
dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)

# Get q vals
qval <- DEXSeq::perGeneQValue(dxr)

# Format q vals
dxr.g <- data.frame(gene = names(qval), qval)

# dxr tsv
keep <- colnames(dxr)[sapply(dxr, class) %in% c("numeric", "character")]
dxr.tsv <- dxr[ , colnames(dxr) %in% keep]

################################
######### Save outputs #########
################################

# dxd
saveRDS(dxd, "dxd_exon.rds")

# results
saveRDS(dxr, "dxr_exon.rds")
write.table(dxr.tsv, "dxr_exon.tsv", sep="\t", quote=FALSE, row.names = TRUE)

# qvals
saveRDS(qval, "qval_exon.rds")
write.table(dxr.g, "dxr_exon.g.tsv", sep="\t", quote=FALSE, row.names = TRUE)

####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("DEXSeq")
sessionInfo()
