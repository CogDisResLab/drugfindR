#' Getters and Setters for drugfindRCoreData slots
#'
#' These functions provide access to the individual slots of a
#' `drugfindRCoreData` object.
#'
#' @param object A `drugfindRCoreData` object.
#' @param value The new value to assign to the slot.
#'
#' @return The value of the slot (for getters) or the modified object (for setters).
#' @name drugfindRCoreData-accessors
#' @keywords internal

# Signature
#' @rdname drugfindRCoreData-accessors
setGeneric("getCoreSignature", function(object) standardGeneric("getCoreSignature"))
#' @rdname drugfindRCoreData-accessors
setMethod("getCoreSignature", "drugfindRCoreData", function(object) object@signature)

#' @rdname drugfindRCoreData-accessors
setGeneric("setCoreSignature<-", function(object, value) standardGeneric("setCoreSignature<-"))
#' @rdname drugfindRCoreData-accessors
setReplaceMethod("setCoreSignature", "drugfindRCoreData", function(object, value) {
    object@signature <- value
    object
})

# Filtered Signature
#' @rdname drugfindRCoreData-accessors
setGeneric("getFilteredSignature", function(object) standardGeneric("getFilteredSignature"))
setMethod("getFilteredSignature", "drugfindRCoreData", function(object) object@filteredSignature)

#' @rdname drugfindRCoreData-accessors
setGeneric("setFilteredSignature<-", function(object, value) standardGeneric("setFilteredSignature<-"))
setReplaceMethod("setFilteredSignature", "drugfindRCoreData", function(object, value) {
    object@filteredSignature <- value
    object
})

# Thresholds
#' @rdname drugfindRCoreData-accessors
setGeneric("getFilterThresholdUp", function(object) standardGeneric("getFilterThresholdUp"))
setMethod("getFilterThresholdUp", "drugfindRCoreData", function(object) object@filterThresholdUp)

#' @rdname drugfindRCoreData-accessors
setGeneric("setFilterThresholdUp<-", function(object, value) standardGeneric("setFilterThresholdUp<-"))
setReplaceMethod("setFilterThresholdUp", "drugfindRCoreData", function(object, value) {
    object@filterThresholdUp <- value
    object
})

#' @rdname drugfindRCoreData-accessors
setGeneric("getFilterThresholdDown", function(object) standardGeneric("getFilterThresholdDown"))
setMethod("getFilterThresholdDown", "drugfindRCoreData", function(object) object@filterThresholdDown)

#' @rdname drugfindRCoreData-accessors
setGeneric("setFilterThresholdDown<-", function(object, value) standardGeneric("setFilterThresholdDown<-"))
setReplaceMethod("setFilterThresholdDown", "drugfindRCoreData", function(object, value) {
    object@filterThresholdDown <- value
    object
})

#' @rdname drugfindRCoreData-accessors
setGeneric("getFilterThresholds", function(object) standardGeneric("getFilterThresholds"))
setMethod("getFilterThresholds", "drugfindRCoreData", function(object) {
    c(object@filterThresholdDown, object@filterThresholdUp)
})

#' @rdname drugfindRCoreData-accessors
setGeneric("setFilterThresholds<-", function(object, value) standardGeneric("setFilterThresholds<-"))
setReplaceMethod("setFilterThresholds", "drugfindRCoreData", function(object, value) {
    if (length(value) == 1) {
        object@filterThresholdDown <- -value
        object@filterThresholdUp <- value
        return(object)
    } else if (length(value) == 2) {
        object@filterThresholdDown <- value[1]
        object@filterThresholdUp <- value[2]
        return(object)
    } else {
        stop("Incorrect number of items in length")
    }
})

# Concordance Limits
#' @rdname drugfindRCoreData-accessors
setGeneric("getConcordanceLimitUp", function(object) standardGeneric("getConcordanceLimitUp"))
setMethod("getConcordanceLimitUp", "drugfindRCoreData", function(object) object@concordanceLimitUp)

#' @rdname drugfindRCoreData-accessors
setGeneric("setConcordanceLimitUp<-", function(object, value) standardGeneric("setConcordanceLimitUp<-"))
setReplaceMethod("setConcordanceLimitUp", "drugfindRCoreData", function(object, value) {
    object@concordanceLimitUp <- value
    object
})

#' @rdname drugfindRCoreData-accessors
setGeneric("getConcordanceLimitDown", function(object) standardGeneric("getConcordanceLimitDown"))
setMethod("getConcordanceLimitDown", "drugfindRCoreData", function(object) object@concordanceLimitDown)

#' @rdname drugfindRCoreData-accessors
setGeneric("setConcordanceLimitDown<-", function(object, value) standardGeneric("setConcordanceLimitDown<-"))
setReplaceMethod("setConcordanceLimitDown", "drugfindRCoreData", function(object, value) {
    object@concordanceLimitDown <- value
    object
})
#' @rdname drugfindRCoreData-accessors
setGeneric("getConcordanceLimits", function(object) standardGeneric("getConcordanceLimits"))
setMethod("getConcordanceLimits", "drugfindRCoreData", function(object) {
    c(object@concordanceLimitDown, object@concordanceLimitUp)
})

#' @rdname drugfindRCoreData-accessors
setGeneric("setConcordanceLimits<-", function(object, value) standardGeneric("setConcordanceLimits<-"))
setReplaceMethod("setConcordanceLimits", "drugfindRCoreData", function(object, value) {
    if (length(value) == 1) {
        object@concordanceLimitDown <- -value
        object@concordanceLimitUp <- value
        return(object)
    } else if (length(value) == 2) {
        object@concordanceLimitDown <- value[1]
        object@concordanceLimitUp <- value[2]
        return(object)
    } else {
        stop("Incorrect number of items in length")
    }
})


# Concordants
#' @rdname drugfindRCoreData-accessors
setGeneric("getUnfilteredConcordants", function(object) standardGeneric("getUnfilteredConcordants"))
setMethod("getUnfilteredConcordants", "drugfindRCoreData", function(object) object@unfilteredConcordants)

#' @rdname drugfindRCoreData-accessors
setGeneric("setUnfilteredConcordants<-", function(object, value) standardGeneric("setUnfilteredConcordants<-"))
setReplaceMethod("setUnfilteredConcordants", "drugfindRCoreData", function(object, value) {
    object@unfilteredConcordants <- value
    object
})

#' @rdname drugfindRCoreData-accessors
setGeneric("getFilteredConcordants", function(object) standardGeneric("getFilteredConcordants"))
setMethod("getFilteredConcordants", "drugfindRCoreData", function(object) object@filteredConcordants)

#' @rdname drugfindRCoreData-accessors
setGeneric("setFilteredConcordants<-", function(object, value) standardGeneric("setFilteredConcordants<-"))
setReplaceMethod("setFilteredConcordants", "drugfindRCoreData", function(object, value) {
    object@filteredConcordants <- value
    object
})
