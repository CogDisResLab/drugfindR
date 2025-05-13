#' @include class-drugfindRCoreData.R
NULL

#' Public Result Result Class
#'
#' The [`drugfindRResult`] S4 class acts as a standardized wrapper for
#' the results generated from a [`drugfindRDataset`] class. This allows
#' easy access to the concordant results and filtering them by specific
#' criteria
#'
#' This class retains the original class of the input signature
#' and converts the result to that input when tested.
#'
#' @slot core An internal `drugfindRCoreData` object
#' @slot source_class A character vector describing the original input class
#'
#' @export
setClass("drugfindRResult",
    slots = c(
        result = "drugfindRCoreData",
        source_class = "character"
    )
)

#' Validity check for the public `drugfindRResult` class
#'
#' This method documents the validity check for the public
#' [`drugfindRResult`] class. This check ensures the created
#' object is not initialized or updated in an invalid manner
#'
#' @details
#' This method returns FALSE in the following conditions:
#'
#' * If the `unfilteredConcordants` slot is not NULL
#' * If the `filteredConcordants` slot is not NULL
#'
#' @param object A [`drugfindRResult`] object.
#'
#' @return A logical value declaring whether the input is valid or not
#'
#' @name validObject.drugfindRResult
setValidity("drugfindRResult", function(object) {
    if (is.null(object@core@unfilteredConcordants)) {
        return("A `drugfindRCoreData` object cannot have the @unfilteredConcordants slot unpopulated")
    } else if (is.null(object@core@filteredConcordants)) {
        return("A `drugfindRCoreData` object cannot have the @filteredConcordants slot unpopulated")
    } else {
        TRUE
    }
})
