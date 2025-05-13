#' @include utilities.R
#' @include class-drugfindRCoreData.R
#' @include drugfindRCoreData.R
NULL

#' External Drug Dataset Class
#'
#' The `drugfindRDataset` S4 class acts as a standardized wrapper for
#' any `data.frame`, `tibble`, or `S4Vectors::DataFrame`, enabling
#' consistent downstream processing in the `drugfindR` package.
#'
#' This class retains the original type of the input so that output
#' can be restored to the same structure after internal transformations.
#'
#' @slot core An internal `drugfindRCoreData` object storing the data
#' as a `tibble`
#' @slot source_class A character vector describing the original input class
#'
#' @export
setClass("drugfindRDataset",
    slots = c(
        core = "drugfindRCoreData",
        source_class = "character"
    ),
    prototype = list(
        core = drugfindRCoreData(tibble::tibble()),
        source_class = "tbl_df"
    )
)

#' Validity check for the public `drugfindRDataset` class
#'
#' This method documents the validity check for the public
#' [`drugfindRDataset`] class. This check ensures the created
#' object is not initialized or updated in an invalid manner
#'
#' @details
#' This method returns FALSE in the following conditions:
#'
#' * If the `unfilteredConcordants` slot is not NULL
#' * If the `filteredConcordants` slot is not NULL
#'
#' @param object A [`drugfindRDataset`] object.
#'
#' @return A logical value declaring whether the input is valid or not
#'
#' @name validObject.drugfindRDataset
setValidity("drugfindRDataset", function(object) {
    if (!is.null(object@core@unfilteredConcordants)) {
        return("A `drugfindRCoreData` object cannot have the @unfilteredConcordants slot populated")
    } else if (!is.null(object@core@filteredConcordants)) {
        return("A `drugfindRCoreData` object cannot have the @filteredConcordants slot populated")
    } else {
        TRUE
    }
})
