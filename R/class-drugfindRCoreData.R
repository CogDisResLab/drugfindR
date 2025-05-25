# Register tbl_df as an S4 class
#' @importClassesFrom tibble tbl_df
setOldClass("tbl_df")

# Define a nullable version for stepwise slot population
setClassUnion("tbl_dfOrNULL", c("tbl_df", "NULL"))

# Define an nullable version for missing values
setClassUnion("characterOrNULL", c("character", "NULL"))

#' Internal Class for Normalized Drug Signature Data
#'
#' This class is used to store a standardized internal representation of
#' drug signature data within the drugfindR pipeline. It includes the original
#' input signature, filtering thresholds, and concordance results. This class
#' is **not** intended for direct user interaction**.
#'
#' @slot coreSignature A tibble containing the original gene expression signature.
#' @slot inputLibrary The library to search the input from.
#' This has to be one of c("OE", "CP", "KD")
#' @slot outputLibrary The library to search the output from.
#' This has to be one of c("OE", "CP", "KD")
#' @slot pairedAnalysis A logical indicating whether a [`pairedAnalysis`]
#' should be performed
#' @slot cellLines A character vector of cell lines to restrict the analysis to
#' @slot filterThresholdUp The upper log fold-change threshold
#' used in filtering.
#' @slot filterThresholdDown The lower log fold-change threshold
#' used in filtering.
#' @slot filteredSignature A tibble of genes after filtering.
#' @slot concordanceLimitUp The upper threshold for compound concordance scores.
#' @slot concordanceLimitDown The lower threshold for compound
#' concordance scores.
#' @slot unfilteredConcordants A tibble of all concordant compounds.
#' @slot filteredConcordants A tibble of filtered concordant compounds.
#'
#' @importFrom tibble tibble
#' @keywords internal
setClass("drugfindRCoreData",
    slots = c(
        coreSignature = "tbl_df",
        inputLibrary = "characterOrNULL",
        outputLibrary = "characterOrNULL",
        pairedAnalysis = "logical",
        cellLines = "characterOrNULL",
        filterThresholdUp = "numeric",
        filterThresholdDown = "numeric",
        filteredSignature = "tbl_dfOrNULL",
        concordanceLimitUp = "numeric",
        concordanceLimitDown = "numeric",
        unfilteredConcordants = "tbl_dfOrNULL",
        filteredConcordants = "tbl_dfOrNULL"
    ),
    prototype = list(
        coreSignature = tibble::tibble(),
        inputLibrary = NULL,
        outputLibrary = NULL,
        pairedAnalysis = TRUE,
        cellLines = NULL,
        filterThresholdUp = 0L,
        filterThresholdDown = 0L,
        filteredSignature = NULL,
        concordanceLimitUp = 0.2,
        concordanceLimitDown = -0.2,
        unfilteredConcordants = NULL,
        filteredConcordants = NULL
    )
)

#' Validity check for the internal `drugfindRCoreData` class
#'
#' This method documents the validity check for the internal
#' [`drugfindRCoreData`] class. This check ensures the created
#' object is not initialized or updated in an invalid manner
#'
#' @details
#' This method returns FALSE in the following conditions:
#'
#' * If the inputLibrary parameters is not one of c("OE", "CP", "KD")
#' * If the outputLibrary parameters is not one of c("OE", "CP", "KD")
#' * If the given cell lines include non-LINCS cell lines
#' * If the concordanceLimits are more than 1 or less than -1
#'
#' @param object A [`drugfindRCoreData`] object.
#'
#' @keywords internal
#' @return A logical value declaring whether the input is valid or not
#'
#' @name validObject.drugfindRCoreData
setValidity("drugfindRCoreData", function(object) {
    if (!validateLibraries(object@inputLibrary)) {
        "@inputLibrary must be one of 'OE', 'KD' or 'CP'"
    } else if (!validateLibraries(object@outputLibrary)) {
        "@outputLibrary must be one of 'OE', 'KD' or 'CP'"
    } else if (!validateCellLines(object@cellLines)) {
        "@cellLines  must be one of the valid cell lines"
    } else if (
        object@concordanceLimitUp > 1L |
            object@concordanceLimitUp < -1L) {
        "@concordanceLimitUp must be between -1 and 1"
    } else if (
        object@concordanceLimitDown > 1L |
            object@concordanceLimitDown < -1L) {
        "@concordanceLimitDown must be between -1 and 1"
    } else {
        TRUE
    }
})
