#!/usr/bin/env Rscript
# Scripts adjusted from Bioconductor and IsoformSwitchAnalyzeR source code
# Please see following for details:
# https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
# and
# https://github.com/kvittingseerup/IsoformSwitchAnalyzeR
# respectivly

# Parse command arguments

argv <- commandArgs(trailingOnly = TRUE)

salmon_output <- getwd()

gtf <- argv[1]

transcript_sequences <- argv[2]

samplesheet <- argv[3]

contrastsheet <- argv[4]

alpha <- as.numeric(argv[5]) # Must be >= 0 and <= 1

dIFcutoff <- as.numeric(argv[6]) # Must be >= 0 and <= 1


# Attach required packages

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))


# Redefine high level functions to add alpha and dIFcutoff

isoformSwitchAnalysisPart1 <- function(

    # Arguments

    switchAnalyzeRlist,
    pathToGTF           = NULL,
    pathToOutput        = NULL,
    alpha               = NULL,
    dIFcutoff           = NULL
) {
    isConditional <- switchAnalyzeRlist$sourceId != 'preDefinedSwitches'

    # preFilter
    if(isConditional) {
        switchAnalyzeRlist <-
            preFilter(
                switchAnalyzeRlist = switchAnalyzeRlist,
                removeSingleIsoformGenes = TRUE,
                quiet = TRUE,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

    }

    # Test isoform switches
    if(isConditional) {

        if(any( switchAnalyzeRlist$conditions$nrReplicates > 5)) {
            switchAnalyzeRlist <-
                isoformSwitchTestSatuRn(
                    switchAnalyzeRlist,
                    reduceToSwitchingGenes = TRUE,
                    quiet = TRUE,
                    alpha = alpha,
                    dIFcutoff = dIFcutoff
                )

        } else {
            switchAnalyzeRlist <-
                isoformSwitchTestDEXSeq(
                    switchAnalyzeRlist,
                    reduceToSwitchingGenes = FALSE,
                    quiet = TRUE,
                    alpha = alpha,
                    dIFcutoff = dIFcutoff
                )
        }

        if (nrow(switchAnalyzeRlist$isoformSwitchAnalysis) == 0) {
            stop('No isoform switches were identified with the current cutoffs.')
        }
    }


    # Predict ORF

    if ( is.null(switchAnalyzeRlist$orfAnalysis) ) {

        # Add known annoation

        switchAnalyzeRlist <- addORFfromGTF(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToGTF = pathToGTF,
            quiet = TRUE
        )

        # Predict novel once (if any are missing)

        if ( any( switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet' )) {
            switchAnalyzeRlist <- analyzeNovelIsoformORF(
                switchAnalyzeRlist = switchAnalyzeRlist,
                analysisAllIsoformsWithoutORF = TRUE,
            )

        }

    }

    # Extract and write sequences

    switchAnalyzeRlist <- extractSequence(
        switchAnalyzeRlist = switchAnalyzeRlist,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        addToSwitchAnalyzeRlist = TRUE,
        writeToFile = FALSE,
        pathToOutput = pathToOutput,
        quiet = TRUE,
        alpha = alpha,
        dIFcutoff = dIFcutoff
    )

    return(switchAnalyzeRlist)
}

isoformSwitchAnalysisPart2 <- function(

    # Core arguments

    switchAnalyzeRlist,

    # Analysis and output arguments

    pathToOutput = NULL,

    # Other arguments

    alpha = NULL,
    dIFcutoff = NULL

) {

    # Predict intron retentions

    switchAnalyzeRlist <-
        analyzeAlternativeSplicing(
            switchAnalyzeRlist = switchAnalyzeRlist,
            quiet = TRUE,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )

    # Predict functional consequences

    switchAnalyzeRlist <-
        analyzeSwitchConsequences(
            switchAnalyzeRlist = switchAnalyzeRlist,
            consequencesToAnalyze = c('intron_retention', 'ORF_seq_similarity', 'NMD_status'),
            quiet = TRUE,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )

    # Make isoform switch plots

    switchPlotTopSwitches(
        switchAnalyzeRlist = switchAnalyzeRlist,
        n = Inf,
        pathToOutput = pathToOutput,
        filterForConsequences = TRUE,
        splitFunctionalConsequences = FALSE,
        quiet = TRUE,
        alpha = alpha,
        dIFcutoff = dIFcutoff
    )

    # Make overall consequences

    pdf(
        file = paste(
            pathToOutput,
            'common_switch_consequences.pdf',
            sep = ''
        ),
        width = 10,
        height = 7
    )
    print(
        extractConsequenceSummary(
            switchAnalyzeRlist = switchAnalyzeRlist,
            plot = TRUE,
            returnResult = FALSE,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )
    )
    dev.off()

    return(switchAnalyzeRlist)
}

isoformSwitchAnalysisCombined <- function(

    # Core arguments

    switchAnalyzeRlist,

    # Annotation arguments

    pathToGTF = NULL,

    # Analysis and output arguments

    pathToOutput = NULL,

    # Other arguments

    alpha = NULL,
    dIFcutoff = NULL
) {

    # Run part 1

    switchAnalyzeRlist <-
        isoformSwitchAnalysisPart1(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToOutput = pathToOutput,
            pathToGTF = pathToGTF,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )

    # Run part 2 without annoation

    switchAnalyzeRlist <-
        isoformSwitchAnalysisPart2(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToOutput = pathToOutput,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )

    return(switchAnalyzeRlist)
}


# Load Salmon output

salmonQuant <- importIsoformExpression(
    parentDir = salmon_output,
    showProgress = FALSE,
    quiet = TRUE
)


# Build desing data frame

design <- read.csv(samplesheet, check.names = FALSE)

design <- design[, c("sample", "condition"), drop = FALSE]

colnames(design) <- c("sampleID", "condition")


# Build contrasts data frame

if (file.exists(contrastsheet)) {

    contrasts <- read.csv(contrastsheet, check.names = FALSE)

    contrasts <- contrasts[, c("treatment", "control"), drop = FALSE]

    colnames(contrasts) <- c("condition_1", "condition_2")

} else {

    contrasts <- NULL

}


# Build swtich list

SwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = design,
    isoformExonAnnoation = gtf,
    isoformNtFasta       = transcript_sequences,
    comparisonsToMake    = contrasts,
    showProgress = FALSE,
    quiet = TRUE
)


tryCatch({

    # Run IsoromSwitchAnalyzR

    SwitchList <- isoformSwitchAnalysisCombined(
        SwitchList,
        pathToGTF = gtf,
        pathToOutput = "results",
        alpha = alpha,
        dIFcutoff = dIFcutoff
    )


    # Save summary

    write.csv(extractSwitchSummary(
        SwitchList,
        filterForConsequences = FALSE
        ),
        "isoformswitchanalyzer_summary.csv"
    )

    # Save isoformFeatures as csv

    write.csv(SwitchList$isoformFeatures, "isoformswitchanalyzer_isoformfeatures.csv")

}, error = function(e) {


    # If an error occurs, write message into summary no index

    write.csv(data.frame(message = conditionMessage(e)),
        "isoformswitchanalyzer_summary.csv",
        row.names = FALSE
    )

    # Write blank isoformFeatures

    write.csv(data.frame(),
        "isoformswitchanalyzer_isoformfeatures.csv",
        row.names = FALSE
    )
})


# Save SwitchList object

saveRDS(SwitchList, "switchlist.rds")


####################################
########### Session info ###########
####################################

# Print sessioninfo to standard out
citation("IsoformSwitchAnalyzeR")
sessionInfo()
