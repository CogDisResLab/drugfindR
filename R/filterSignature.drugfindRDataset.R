#' @include specialGenerics.R
NULL

#' Filter the signature to the defined parameters
#'
#' This function performs the filtering step. The filtering
#' is based on the logFC threshold set on initialization.
#'
#' The default filtering value is 0, which means that no
#' filtering is done.
#'
#' @param object An object of type [`drugfindRDataset`]
#'
#' @return A copy of the input [`drugfindRDataset`], with the @filteredSignature
#' slot populated
#'
#' @importFrom methods validObject
#' @importFrom dplyr filter
#'
#' @export
setMethod("filterSignature", "drugfindRDataset", function(object) {
    filteredObject <- filterSignature(object@core)
    object@core <- filteredobject
    validObject(object)
    object
})
