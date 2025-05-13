#' @include class-drugfindRCoreData.R
NULL

#' Accessor Methods for drugfindRCoreData Slots
#'
#' These functions provide internal access to the individual slots of a
#' `drugfindRCoreData` object used in the drugfindR pipeline.
#'
#' @param object A [`drugfindRCoreData`] object.
#' @param value The value to assign to a slot (for setters only).
#'
#' @return For getters, the value stored in the relevant slot. For setters, the updated object.
#'
#' ## Getter methods
#'
#' ### Input Parameters
#' - `getCoreSignature(object)` returns the core gene signature.
#' - `getFilterThresholds(object)` returns both thresholds as a vector.
#' - `getConcordanceLimits(object)` returns both concordance limits as a vector.
#'
#' ### Intermediate and Final Results
#' - `getFilteredSignature(object)` returns the filtered signature.
#' - `getUnfilteredConcordants(object)` returns unfiltered concordant drugs.
#' - `getFilteredConcordants(object)` returns filtered concordant drugs.
#'
#' ## Setter methods
#'
#' ### Input Parameters
#' - `setFilterThresholdUp(object) <- value` sets the upper threshold.
#' - `setFilterThresholdDown(object) <- value` sets the lower threshold.
#' - `setFilterThresholds(object) <- value` sets both thresholds from a vector of length 1 or 2.
#' - `setConcordanceLimitUp(object) <- value` sets the upper limit.
#' - `setConcordanceLimitDown(object) <- value` sets the lower limit.
#' - `setConcordanceLimits(object) <- value` sets both concordance limits from a vector of length 1 or 2.
#' ### Intermediate and Final Results
#' - `setFilteredSignature(object) <- value` sets the filtered signature.
#' - `setFilteredConcordants(object) <- value` sets filtered concordants.
#'
#' @name drugfindRCoreData-accessors
#' @aliases getCoreSignature getFilteredSignature getFilterThresholdUp
#' getFilterThresholdDown getFilterThresholds getConcordanceLimitUp
#' getConcordanceLimitDown getConcordanceLimits
#' setFilterThresholdUp<- setFilterThresholdDown<- setFilterThresholds<-
#' setConcordanceLimitUp<- setConcordanceLimitDown<- setConcordanceLimits<-
#' getUnfilteredConcordants getFilteredConcordants
#'
#' @importFrom methods validObject
#'
#' @keywords internal

### -- Signature --
setGeneric("getCoreSignature", function(object) standardGeneric("getCoreSignature"))
setMethod("getCoreSignature", "drugfindRCoreData", function(object) object@signature)

### -- Filtered Signature --
setGeneric("getFilteredSignature", function(object) standardGeneric("getFilteredSignature"))
setMethod("getFilteredSignature", "drugfindRCoreData", function(object) object@filteredSignature)

### -- Filter Thresholds --
setGeneric("getFilterThresholdUp", function(object) standardGeneric("getFilterThresholdUp"))
setMethod("getFilterThresholdUp", "drugfindRCoreData", function(object) object@filterThresholdUp)

setGeneric("setFilterThresholdUp<-", function(object, value) standardGeneric("setFilterThresholdUp<-"))
setReplaceMethod("setFilterThresholdUp", "drugfindRCoreData", function(object, value) {
    object@filterThresholdUp <- value
    validObject(object)
    object
})

setGeneric("getFilterThresholdDown", function(object) standardGeneric("getFilterThresholdDown"))
setMethod("getFilterThresholdDown", "drugfindRCoreData", function(object) object@filterThresholdDown)

setGeneric("setFilterThresholdDown<-", function(object, value) standardGeneric("setFilterThresholdDown<-"))
setReplaceMethod("setFilterThresholdDown", "drugfindRCoreData", function(object, value) {
    object@filterThresholdDown <- value
    validObject(object)
    object
})

setGeneric("getFilterThresholds", function(object) standardGeneric("getFilterThresholds"))
setMethod("getFilterThresholds", "drugfindRCoreData", function(object) {
    c(getFilterThresholdDown(object), getFilterThresholdUp(object))
})

setGeneric("setFilterThresholds<-", function(object, value) standardGeneric("setFilterThresholds<-"))
setReplaceMethod("setFilterThresholds", "drugfindRCoreData", function(object, value) {
    if (length(value) == 1) {
        object@filterThresholdDown <- -value
        object@filterThresholdUp <- value
    } else if (length(value) == 2) {
        object@filterThresholdDown <- value[1]
        object@filterThresholdUp <- value[2]
    } else {
        stop("Incorrect number of items in length")
    }
    validObject(object)
    object
})

### -- Concordance Limits --
setGeneric("getConcordanceLimitUp", function(object) standardGeneric("getConcordanceLimitUp"))
setMethod("getConcordanceLimitUp", "drugfindRCoreData", function(object) object@concordanceLimitUp)

setGeneric("setConcordanceLimitUp<-", function(object, value) standardGeneric("setConcordanceLimitUp<-"))
setReplaceMethod("setConcordanceLimitUp", "drugfindRCoreData", function(object, value) {
    object@concordanceLimitUp <- value
    validObject(object)
    object
})

setGeneric("getConcordanceLimitDown", function(object) standardGeneric("getConcordanceLimitDown"))
setMethod("getConcordanceLimitDown", "drugfindRCoreData", function(object) object@concordanceLimitDown)

setGeneric("setConcordanceLimitDown<-", function(object, value) standardGeneric("setConcordanceLimitDown<-"))
setReplaceMethod("setConcordanceLimitDown", "drugfindRCoreData", function(object, value) {
    object@concordanceLimitDown <- value
    validObject(object)
    object
})

setGeneric("getConcordanceLimits", function(object) standardGeneric("getConcordanceLimits"))
setMethod("getConcordanceLimits", "drugfindRCoreData", function(object) {
    c(getConcordanceLimitDown(object), getConcordanceLimitUp(object))
})

setGeneric("setConcordanceLimits<-", function(object, value) standardGeneric("setConcordanceLimits<-"))
setReplaceMethod("setConcordanceLimits", "drugfindRCoreData", function(object, value) {
    if (length(value) == 1) {
        object@concordanceLimitDown <- -value
        object@concordanceLimitUp <- value
    } else if (length(value) == 2) {
        object@concordanceLimitDown <- value[1]
        object@concordanceLimitUp <- value[2]
    } else {
        stop("Incorrect number of items in length")
    }
    validObject(object)
    object
})

### -- Concordants --
setGeneric("getUnfilteredConcordants", function(object) standardGeneric("getUnfilteredConcordants"))
setMethod("getUnfilteredConcordants", "drugfindRCoreData", function(object) object@unfilteredConcordants)

setGeneric("getFilteredConcordants", function(object) standardGeneric("getFilteredConcordants"))
setMethod("getFilteredConcordants", "drugfindRCoreData", function(object) object@filteredConcordants)
