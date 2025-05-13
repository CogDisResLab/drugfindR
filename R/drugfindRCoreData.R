#' @include class-drugfindRCoreData.R
NULL

#' Construct a drugfindRCoreData Object
#'
#' This internal constructor creates a standardized signature container
#' for use in the drugfindR pipeline. Most fields are optional and may be filled
#' incrementally during filtering and concordance analysis.
#'
#' @param signature A tibble containing the original gene expression data.
#' @param inputLibrary The library to search the input from. One of c("OE", "CP", "KD")
#' @param outputLibrary The library to search the output from. One of c("OE", "CP", "KD")
#' @param pairedAnalysis A logical indicating whether a [`pairedAnalysis`] should be performed
#' @param cellLines A character vector of cell lines to restrict the analysis to
#' @param filterThresholdUp Upper logFC threshold used in filtering.
#' @param filterThresholdDown Lower logFC threshold used in filtering.
#' @param concordanceLimitUp Upper threshold for selecting concordant compounds.
#' @param concordanceLimitDown Lower threshold for selecting concordant compounds.
#'
#' @return A `drugfindRCoreData` S4 object
#' @keywords internal
drugfindRCoreData <- function(signature,
                              inputLibrary = NULL,
                              outputLibrary = NULL,
                              pairedAnalysis = TRUE,
                              cellLines = NULL,
                              filterThresholdUp = 0,
                              filterThresholdDown = 0,
                              concordanceLimitUp = 0.2,
                              concordanceLimitDown = -0.2) {
    methods::new("drugfindRCoreData",
        signature = signature,
        inputLibrary = inputLibrary,
        outputLibrary = outputLibrary,
        pairedAnalysis = pairedAnalysis,
        cellLines = cellLines,
        filterThresholdUp = filterThresholdUp,
        filterThresholdDown = filterThresholdDown,
        filteredSignature = NULL,
        concordanceLimitUp = concordanceLimitUp,
        concordanceLimitDown = concordanceLimitDown,
        unfilteredConcordants = NULL,
        filteredConcordants = NULL
    )
}
