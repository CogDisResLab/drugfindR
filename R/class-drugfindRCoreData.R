# Register tbl_df as an S4 class
#' @importClassesFrom tibble tbl_df
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
