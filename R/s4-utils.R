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
setMethod("dimnames", "drugfindRDataset", function(x) {
    dimnames(x@core)
})

# #' Head for a drugfindRDataset
# #'
# #' @param x A drugfindRDataset object
# #' @param ... Additional arguments passed to utils::head
# #' @return The first rows of the signature data as tibble
# #' @export
# setMethod("head", "drugfindRDataset", function(x, ...) {
#     utils::head(x@signature@core@data, ...)
# })
#
# #' Number of rows for a drugfindRResult
# #'
# #' @param x A drugfindRResult object
# #' @return Integer number of rows in the filtered concordants
# #' @export
# setMethod("nrow", "drugfindRResult", function(x) {
#     nrow(x@filteredConcordants)
# })
#
# #' Number of columns for a drugfindRResult
# #'
# #' @param x A drugfindRResult object
# #' @return Integer number of columns in the filtered concordants
# #' @export
# setMethod("ncol", "drugfindRResult", function(x) {
#     ncol(x@filteredConcordants)
# })
#
# #' Dimensions for a drugfindRResult
# #'
# #' @param x A drugfindRResult object
# #' @return Integer vector: rows and columns of filtered concordants
# #' @export
# setMethod("dim", "drugfindRResult", function(x) {
#     dim(x@filteredConcordants)
# })
#
# #' Head for a drugfindRResult
# #'
# #' @param x A drugfindRResult object
# #' @param ... Additional args passed to utils::head
# #' @return The first rows of filtered concordants
# #' @export
# setMethod("head", "drugfindRResult", function(x, ...) {
#     utils::head(x@filteredConcordants, ...)
# })
