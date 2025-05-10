#' @include class-drugfindRCoreData.R
NULL

# Preparing to implement `summary` for drugfindRCoreData
setGeneric("summary",
    def = function(object) standardGeneric("summary")
)

#' Summarize a drugfindRCoreData Object
#'
#' This function returns a named list containing all slots from a
#' `drugfindRCoreData` object. It is useful for debugging or inspection.
#'
#' @param object A `drugfindRCoreData` object.
#'
#' @return A named list with each slot's current value.
#' @export
setMethod("summary", "drugfindRCoreData", function(object) {
    stopifnot(inherits(object, "drugfindRCoreData"))
    list(
        signature = object@signature,
        filterThresholdUp = object@filterThresholdUp,
        filterThresholdDown = object@filterThresholdDown,
        filteredSignature = object@filteredSignature,
        concordanceLimitUp = object@concordanceLimitUp,
        concordanceLimitDown = object@concordanceLimitDown,
        unfilteredConcordants = object@unfilteredConcordants,
        filteredConcordants = object@filteredConcordants
    )
})
