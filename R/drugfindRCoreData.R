#' @include class-drugfindRCoreData.R
NULL

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
