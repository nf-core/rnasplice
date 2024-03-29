#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics
# Email: bioinformatics@zifornd.com
# License: MIT

# Splits any input file (e.g. tpm) by column using samplesheet information
# Utility script for Suppa processing
# Also takes note of how many samples for each condition for down stream Clustering module

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {

    stop("Usage: suppa_split_file.R <input_file> <samplesheet> <output_file_suffix> <calculate_ranges> <prefix>", call.=FALSE)

}

######################################
########### Collect inputs ###########
######################################

input_file         <- args[1]
samplesheet        <- args[2]
output_file_suffix <- args[3]
calculate_ranges   <- args[4] # TRUE if we want to return cluster ranges

if (length(args) == 5){

    prefix <- args[5]

} else {

    prefix <- ""

}

######################################
######### Read in input file #########
######################################

input_data <- read.csv(input_file, sep = "\t", header = TRUE, check.names = FALSE)

######################################
####### Process samplesheet ##########
######################################

# Read in samplesheet
samplesheet <- read.csv(samplesheet, header = TRUE, check.names = FALSE)

# check header of sample sheet
if (!c("sample") %in% colnames(samplesheet) | !c("condition") %in% colnames(samplesheet)) {

    stop("suppa_split_file.R Samplesheet must contain 'sample' and 'condition' column headers.", call.=FALSE)

}

# Take only sample and condition columns
samplesheet <- samplesheet[,c("sample", "condition")]

# filter for unique rows based on sample name
samplesheet <- samplesheet[!duplicated(samplesheet[,"sample"]),]

# Take unique conditions
conditions <- unique(samplesheet[,"condition"])

#########################################################
####### Define function for splitting input files #######
#########################################################

# Function for taking all sample names associated with a given condition
split_files <- function(condition, samplesheet, input_data, output_file_suffix, prefix, calculate_ranges){

    # Get indices of rows which cover given condition for ranges
    indices <- which(samplesheet$condition == condition)

    # Get sample names for given condition
    sample_names <- samplesheet[samplesheet$condition == condition,]$sample

    # Check header of input_data contains all samples from processed samplesheet
    if (!all(samplesheet$sample %in% colnames(input_data))) {

        stop("suppa_split_file.R Input_file must contain samplesheet samples.", call.=FALSE)

    }

    if (prefix == ""){

        output_file <- paste0(condition, output_file_suffix)

    } else {

        output_file <- paste0(prefix, "_" ,condition, output_file_suffix)

    }

    # Subset input files and save out as new file
    write.table(input_data[,sample_names, drop=F], file = output_file, quote = FALSE, sep = "\t")

    # Get Cluster ranges which match the tpm and psi files above (1-3 4-6)
    # Column numbers have to be continuous, with no overlapping or missing columns between them. Ex:1-3,4-6
    if (calculate_ranges) {

        if (indices[1] == 1) {

            range <- paste0(as.character(indices[1]), "-" , as.character(indices[length(indices)]))
            cat(range, file = "ranges.txt", sep = "")

        } else {

            range <- paste0(as.character(indices[1]), "-", as.character(indices[length(indices)]))
            cat(c(",",range), file = "ranges.txt", sep = "", append=TRUE)
        }
    }
}

####################################################################
##### Iterate through conditions - split files and get ranges ######
####################################################################

for (cond in conditions) {

    # Split files
    split_files(cond, samplesheet, input_data, output_file_suffix, prefix, calculate_ranges)

}

####################################
########## Session info ############
####################################

# Print sessioninfo to standard out
sessionInfo()

