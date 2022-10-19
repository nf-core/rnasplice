#!/usr/bin/env Rscript

# Splits any input file (e.g. tpm) by column using samplesheet information
# Utility script for Suppa processing 
# Also takes note of how many samples for each condition for down stream Clustering module

# Parse command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  
  stop("Usage: suppa_split_file.R <input_file> <samplesheet> <output_file_suffix> <calculate_ranges>", call.=FALSE)
  
} 

######################################
########### Collect inputs ###########
######################################

input_file <- args[1]
samplesheet <- args[2]
output_file_suffix <- args[3]
calculate_ranges <- args[4] # true if we want to return cluster ranges

######################################
####### Process samplesheet ##########
######################################

# Read in samplesheet
samplesheet <- read.csv(samplesheet, header = TRUE)

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

# Function for taking all sample names associated with a given condition
get_sample_names <- function(condition, samplesheet){

  sample_names <- samplesheet[samplesheet$condition == condition,]$sample
  return(sample_names)
}

# Loop over all unique conditions and retrieve samples for each condition
# set as data frame
samples_cond <- as.data.frame(sapply(conditions, get_sample_names, samplesheet = samplesheet))

######################################
####### Define output files ##########
######################################

# Ready output file names based on user suffix input and unique conditions from samplesheet
output_files <- sapply(colnames(samples_cond), function(x, suffix){paste0(x, suffix)}, suffix = output_file_suffix)

################################################################
########## Process input and save new output files #############
################################################################

# Read in input file
input_file <- read.csv(input_file, sep="\t", header=TRUE)

# Check header of input_file contains all samples from processed samplesheet
if (!all(samplesheet$sample %in% colnames(input_file))) {
  
  stop("suppa_split_file.R Input_file must contain samplesheet samples.", call.=FALSE)
  
}

# Loop through conditions and create separate output files per unique condition
# Also take note of sample number whilst we are here for down stream clustering module
idx <- 0
range <- ""

for (cond in conditions) {
  
  # Write output files per condition (e.g. tpm and psi files)
  sample_names <- samples_cond[,cond]
  output_file_name <- as.character(output_files[cond])
  write.table(input_file[,sample_names], file = output_file_name, quote = FALSE, sep = "\t")
  
  if (calculate_ranges){

    # Get Cluster ranges which match the tpm and psi files above (1-3 4-6)
    if (idx == 0) {

        range <- paste0("1-" , as.character(length(sample_names)))
        cat(range, file = "ranges.txt", sep = "")
        idx <- 1

    } else {

        prior_group_sum <- as.numeric(substr(range, nchar(range), nchar(range)))
        range <- paste0(as.character(prior_group_sum + 1), "-", as.character(length(sample_names) + prior_group_sum))
        cat(c(",",range), file = "ranges.txt", sep = "", append=TRUE)

    }
  }

}

####################################
########## Session info ############
####################################

# Print sessioninfo to standard out
sessionInfo()


