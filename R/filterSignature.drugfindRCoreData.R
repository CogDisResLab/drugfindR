#' @include utilities.R
NULL

### Method to populate filtered signature slot
setGeneric(
    "filterSignature",
    function(object) standardGeneric("filterSignature")
)
#' Filter the signature to the defined parameters
#'
#' This function performs the filtering step. The filtering
#' is based on the logFC threshold set on initialization.
#'
#' The default filtering value is 0, which means that no
#' filtering is done.
#'
#' @param object An object of type [`drugfindRCoreData`]
#'
#' @return A copy of the input object, with the @filteredSignature
#' slot populated with the filteredSignature
#'
#' @importFrom methods validObject
#' @importFrom dplyr filter
#'
#' @export
setMethod("filterSignature", "drugfindRCoreData", function(object) {
    filtered <- object@signature |>
        dplyr::filter(
            Value_LogDiffExp > getFilterThresholdUp(object) |
                Value_LogDiffExp < getFilterThresholdDown(object)
        )
    object@filteredSignature <- filtered
    validObject(object)
    object
})
