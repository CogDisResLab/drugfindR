# Register tbl_df as an S4 class
setOldClass("tbl_df")

# Define a nullable version for stepwise slot population
setClassUnion("tbl_dfOrNULL", c("tbl_df", "NULL"))

#' Internal Class for Normalized Drug Signature Data
#'
#' This class is used to store a standardized internal representation of
#' drug signature data within the drugfindR pipeline. It includes the original
#' input signature, filtering thresholds, and concordance results. This class
#' is **not intended for direct user interaction**.
#'
#' @slot signature A tibble containing the original gene expression signature.
#' @slot filterThresholdUp The upper log fold-change threshold used in filtering.
#' @slot filterThresholdDown The lower log fold-change threshold used in filtering.
#' @slot filteredSignature A tibble of genes after filtering.
#' @slot concordanceLimitUp The upper threshold for compound concordance scores.
#' @slot concordanceLimitDown The lower threshold for compound concordance scores.
#' @slot unfilteredConcordants A tibble of all concordant compounds.
#' @slot filteredConcordants A tibble of filtered concordant compounds.
#'
#' @importFrom tibble tibble
#' @keywords internal
setClass("drugfindRCoreData",
    slots = c(
        signature = "tbl_df",
        filterThresholdUp = "numeric",
        filterThresholdDown = "numeric",
        filteredSignature = "tbl_dfOrNULL",
        concordanceLimitUp = "numeric",
        concordanceLimitDown = "numeric",
        unfilteredConcordants = "tbl_dfOrNULL",
        filteredConcordants = "tbl_dfOrNULL"
    ),
    prototype = list(
        signature = tibble::tibble(),
        filterThresholdUp = NA_real_,
        filterThresholdDown = NA_real_,
        filteredSignature = NULL,
        concordanceLimitUp = NA_real_,
        concordanceLimitDown = NA_real_,
        unfilteredConcordants = NULL,
        filteredConcordants = NULL
    )
)

#' Construct a drugfindRCoreData Object
#'
#' This internal constructor creates a standardized signature container
#' for use in the drugfindR pipeline. Most fields are optional and may be filled
#' incrementally during filtering and concordance analysis.
#'
#' @param signature A tibble containing the original gene expression data.
#' @param filterThresholdUp Upper logFC threshold used in filtering.
#' @param filterThresholdDown Lower logFC threshold used in filtering.
#' @param filteredSignature A filtered tibble of genes meeting the threshold.
#' @param concordanceLimitUp Upper threshold for selecting concordant compounds.
#' @param concordanceLimitDown Lower threshold for selecting concordant compounds.
#' @param unfilteredConcordants All raw compound concordance scores (optional).
#' @param filteredConcordants Filtered set of top compound concordants (optional).
#'
#' @return A `drugfindRCoreData` S4 object
#' @keywords internal
#'
#' @examples
#' sig <- tibble::tibble(Gene = c("A", "B"), Value_LogDiffExp = c(-1, 1))
#' core <- drugfindRCoreData(signature = sig)
drugfindRCoreData <- function(signature,
                              filterThresholdUp = NA_real_,
                              filterThresholdDown = NA_real_,
                              filteredSignature = NULL,
                              concordanceLimitUp = NA_real_,
                              concordanceLimitDown = NA_real_,
                              unfilteredConcordants = NULL,
                              filteredConcordants = NULL) {
    if (!inherits(signature, "tbl_df")) {
        stop("`signature` must be a tibble.")
    }

    methods::new("drugfindRCoreData",
        signature = signature,
        filterThresholdUp = filterThresholdUp,
        filterThresholdDown = filterThresholdDown,
        filteredSignature = filteredSignature,
        concordanceLimitUp = concordanceLimitUp,
        concordanceLimitDown = concordanceLimitDown,
        unfilteredConcordants = unfilteredConcordants,
        filteredConcordants = filteredConcordants
    )
}
