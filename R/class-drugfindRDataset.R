#' @include class-drugfindRCoreData.R
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
    )
)
