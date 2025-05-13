#' @include class-drugfindRDataset.R
#' @include class-drugfindRCoreData.R
#' @include class-drugfindRResult.R
NULL

#' Dimensions for a drugfindRCoreData
#'
#' @param x A drugfindRCoreData x
#' @return Integer vector: rows and columns of the signature data
#' @export
setMethod("dim", "drugfindRCoreData", function(x) {
    dim(x@signature)
})

#' Dimension names for a drugfindRCoreData
#' @param x A drugfindRCoreData x
#' @return Character vector: rows and columns of the signature data
#' @export
setMethod("dimnames", "drugfindRCoreData", function(x) {
    dimnames(x@signature)
})

#' Column names for a drugfindRCoreData
#' @param x A drugfindRCoreData x
#' @return Character vector: rows and columns of the signature data
#' @export
setMethod("names", "drugfindRCoreData", function(x) {
    names(x@signature)
})

#' Dimensions for a drugfindRDataset
#'
#' @param x A drugfindRDataset x
#' @return Integer vector: rows and columns of the signature data
#' @export
setMethod("dim", "drugfindRDataset", function(x) {
    dim(x@core)
})

#' Dimension names for a drugfindRDataset
#'
#' @param x A drugfindRDataset x
#' @return Character vector: rows and columns of the signature data
#' @export
setMethod("dimnames", "drugfindRDataset", function(x) {
    dimnames(x@core)
})

#' Dimension names for a drugfindRDataset
#'
#' @param x A drugfindRDataset x
#' @return Character vector: rows and columns of the signature data
#' @export
setMethod("dim", "drugfindRDataset", function(x) {
    dim(x@core)
})

#' Number of columns for a drugfindRResult object
#'
#' @param x A drugfindRResult object
#' @return Integer number of columns in the filtered concordants
#' @export
setMethod("dimnames", "drugfindRResult", function(x) {
    dimnames(x@filteredConcordants)
})

#' Dimensions for a drugfindRResult
#'
#' @param x A drugfindRResult object
#' @return Integer vector: rows and columns of filtered concordants
#' @export
setMethod("names", "drugfindRResult", function(x) {
    names(x@filteredConcordants)
})
