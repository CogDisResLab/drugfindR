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
#' @section Core Signature:
#' - `getCoreSignature(object)` returns the core gene signature.
#' - `setCoreSignature(object) <- value` sets the core gene signature.
#'
#' @section Filtered Signature:
#' - `getFilteredSignature(object)` returns the filtered signature.
#' - `setFilteredSignature(object) <- value` sets the filtered signature.
#'
#' @section Filter Thresholds:
#' - `getFilterThresholdUp(object)` returns the upper filter threshold.
#' - `setFilterThresholdUp(object) <- value` sets the upper threshold.
#' - `getFilterThresholdDown(object)` returns the lower filter threshold.
#' - `setFilterThresholdDown(object) <- value` sets the lower threshold.
#' - `getFilterThresholds(object)` returns both thresholds as a vector.
#' - `setFilterThresholds(object) <- value` sets both thresholds from a vector of length 1 or 2.
#'
#' @section Concordance Limits:
#' - `getConcordanceLimitUp(object)` returns the upper concordance limit.
#' - `setConcordanceLimitUp(object) <- value` sets the upper limit.
#' - `getConcordanceLimitDown(object)` returns the lower concordance limit.
#' - `setConcordanceLimitDown(object) <- value` sets the lower limit.
#' - `getConcordanceLimits(object)` returns both concordance limits as a vector.
#' - `setConcordanceLimits(object) <- value` sets both concordance limits from a vector of length 1 or 2.
#'
#' @section Concordants:
#' - `getUnfilteredConcordants(object)` returns unfiltered concordant drugs.
#' - `setUnfilteredConcordants(object) <- value` sets unfiltered concordants.
#' - `getFilteredConcordants(object)` returns filtered concordant drugs.
#' - `setFilteredConcordants(object) <- value` sets filtered concordants.
#'
#' @name drugfindRCoreData-accessors
#' @aliases getCoreSignature setCoreSignature<- getFilteredSignature setFilteredSignature<-
#'   getFilterThresholdUp setFilterThresholdUp<- getFilterThresholdDown setFilterThresholdDown<-
#'   getFilterThresholds setFilterThresholds<- getConcordanceLimitUp setConcordanceLimitUp<-
#'   getConcordanceLimitDown setConcordanceLimitDown<- getConcordanceLimits setConcordanceLimits<-
#'   getUnfilteredConcordants setUnfilteredConcordants<- getFilteredConcordants setFilteredConcordants<-
#' @keywords internal

### -- Signature --
setGeneric("getCoreSignature", function(object) standardGeneric("getCoreSignature"))
setMethod("getCoreSignature", "drugfindRCoreData", function(object) object@signature)

setGeneric("setCoreSignature<-", function(object, value) standardGeneric("setCoreSignature<-"))
setReplaceMethod("setCoreSignature", "drugfindRCoreData", function(object, value) {
    object@signature <- value
    object
})

### -- Filtered Signature --
setGeneric("getFilteredSignature", function(object) standardGeneric("getFilteredSignature"))
setMethod("getFilteredSignature", "drugfindRCoreData", function(object) object@filteredSignature)

setGeneric("setFilteredSignature<-", function(object, value) standardGeneric("setFilteredSignature<-"))
setReplaceMethod("setFilteredSignature", "drugfindRCoreData", function(object, value) {
    object@filteredSignature <- value
    object
})

### -- Filter Thresholds --
setGeneric("getFilterThresholdUp", function(object) standardGeneric("getFilterThresholdUp"))
setMethod("getFilterThresholdUp", "drugfindRCoreData", function(object) object@filterThresholdUp)

setGeneric("setFilterThresholdUp<-", function(object, value) standardGeneric("setFilterThresholdUp<-"))
setReplaceMethod("setFilterThresholdUp", "drugfindRCoreData", function(object, value) {
    object@filterThresholdUp <- value
    object
})

setGeneric("getFilterThresholdDown", function(object) standardGeneric("getFilterThresholdDown"))
setMethod("getFilterThresholdDown", "drugfindRCoreData", function(object) object@filterThresholdDown)

setGeneric("setFilterThresholdDown<-", function(object, value) standardGeneric("setFilterThresholdDown<-"))
setReplaceMethod("setFilterThresholdDown", "drugfindRCoreData", function(object, value) {
    object@filterThresholdDown <- value
    object
})

setGeneric("getFilterThresholds", function(object) standardGeneric("getFilterThresholds"))
setMethod("getFilterThresholds", "drugfindRCoreData", function(object) {
    c(object@filterThresholdDown, object@filterThresholdUp)
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
    object
})

### -- Concordance Limits --
setGeneric("getConcordanceLimitUp", function(object) standardGeneric("getConcordanceLimitUp"))
setMethod("getConcordanceLimitUp", "drugfindRCoreData", function(object) object@concordanceLimitUp)

setGeneric("setConcordanceLimitUp<-", function(object, value) standardGeneric("setConcordanceLimitUp<-"))
setReplaceMethod("setConcordanceLimitUp", "drugfindRCoreData", function(object, value) {
    object@concordanceLimitUp <- value
    object
})

setGeneric("getConcordanceLimitDown", function(object) standardGeneric("getConcordanceLimitDown"))
setMethod("getConcordanceLimitDown", "drugfindRCoreData", function(object) object@concordanceLimitDown)

setGeneric("setConcordanceLimitDown<-", function(object, value) standardGeneric("setConcordanceLimitDown<-"))
setReplaceMethod("setConcordanceLimitDown", "drugfindRCoreData", function(object, value) {
    object@concordanceLimitDown <- value
    object
})

setGeneric("getConcordanceLimits", function(object) standardGeneric("getConcordanceLimits"))
setMethod("getConcordanceLimits", "drugfindRCoreData", function(object) {
    c(object@concordanceLimitDown, object@concordanceLimitUp)
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
    object
})

### -- Concordants --
setGeneric("getUnfilteredConcordants", function(object) standardGeneric("getUnfilteredConcordants"))
setMethod("getUnfilteredConcordants", "drugfindRCoreData", function(object) object@unfilteredConcordants)

setGeneric("setUnfilteredConcordants<-", function(object, value) standardGeneric("setUnfilteredConcordants<-"))
setReplaceMethod("setUnfilteredConcordants", "drugfindRCoreData", function(object, value) {
    object@unfilteredConcordants <- value
    object
})

setGeneric("getFilteredConcordants", function(object) standardGeneric("getFilteredConcordants"))
setMethod("getFilteredConcordants", "drugfindRCoreData", function(object) object@filteredConcordants)

setGeneric("setFilteredConcordants<-", function(object, value) standardGeneric("setFilteredConcordants<-"))
setReplaceMethod("setFilteredConcordants", "drugfindRCoreData", function(object, value) {
    object@filteredConcordants <- value
    object
})
