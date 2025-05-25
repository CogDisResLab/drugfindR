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
#' @return For getters, the value stored in the relevant slot.
#' For setters, the updated object.
#'
#' ## Getter methods
#'
#' ### Input Parameters
#' - `coreSignature(object)` returns the core gene signature.
#' - `filterThresholds(object)` returns both thresholds as a vector.
#' - `concordanceLimits(object)` returns both concordance limits as a vector.
#'
#' ### Intermediate and Final Results
#' - `filteredSignature(object)` returns the filtered signature.
#' - `unfilteredConcordants(object)` returns unfiltered concordant drugs.
#' - `filteredConcordants(object)` returns filtered concordant drugs.
#'
#' ## Setter methods
#'
#' ### Input Parameters
#' - `filterThresholdUp(object) <- value` sets the upper threshold.
#' - `filterThresholdDown(object) <- value` sets the lower threshold.
#' - `filterThresholds(object) <- value` sets both thresholds from a vector
#' of length 1 or 2.
#' - `concordanceLimitUp(object) <- value` sets the upper limit.
#' - `concordanceLimitDown(object) <- value` sets the lower limit.
#' - `concordanceLimits(object) <- value` sets both concordance limits
#' from a vector of length 1 or 2.
#' ### Intermediate and Final Results
#' - `filteredSignature(object) <- value` sets the filtered signature.
#' - `filteredConcordants(object) <- value` sets filtered concordants.
#'
#' @name drugfindRCoreData-accessors
#' @aliases coreSignature filteredSignature filterThresholdUp
#' filterThresholdDown filterThresholds concordanceLimitUp
#' concordanceLimitDown concordanceLimits
#' filterThresholdUp<- filterThresholdDown<- filterThresholds<-
#' concordanceLimitUp<- concordanceLimitDown<- concordanceLimits<-
#' unfilteredConcordants filteredConcordants
#'
#' @importFrom methods validObject
#'
#' @keywords internal

### -- Core Signature --

#### Generics

setGeneric(
    "coreSignature",
    function(object, ...) standardGeneric("coreSignature")
)

setGeneric(
    "setCoreSignature",
    function(object, value) standardGeneric("setCoreSignature")
)

setGeneric(
    "coreSignature<-",
    function(object, value) standardGeneric("coreSignature<-")
)


#### Getter
setMethod(
    "coreSignature",
    "drugfindRCoreData",
    function(object) object@coreSignature
)

#### Setter

setMethod(
    "setCoreSignature",
    "drugfindRCoreData",
    function(object, value) {
        coreSignature(object) <- value
    }
)

#### Replacer

setMethod(
    "coreSignature<-",
    "drugfindRCoreData",
    function(object, value) {
        object@coreSignature <- value
        validObject(object)
        object
    }
)


### -- Input Library --

#### Generics

setGeneric(
    "inputLibrary",
    function(object, ...) standardGeneric("inputLibrary")
)

setGeneric(
    "setInputLibrary",
    function(object, value) standardGeneric("setInputLibrary")
)

setGeneric(
    "inputLibrary<-",
    function(object, value) standardGeneric("inputLibrary<-")
)


#### Getter
setMethod(
    "inputLibrary",
    "drugfindRCoreData",
    function(object) object@inputLibrary
)

#### Setter

setMethod(
    "setInputLibrary",
    "drugfindRCoreData",
    function(object, value) {
        inputLibrary(object) <- value
    }
)

#### Replacer

setMethod(
    "inputLibrary<-",
    "drugfindRCoreData",
    function(object, value) {
        object@inputLibrary <- value
        validObject(object)
        object
    }
)

### -- output Library --

#### Generics

setGeneric(
    "outputLibrary",
    function(object, ...) standardGeneric("outputLibrary")
)

setGeneric(
    "setOutputLibrary",
    function(object, value) standardGeneric("setOutputLibrary")
)

setGeneric(
    "outputLibrary<-",
    function(object, value) standardGeneric("outputLibrary<-")
)


#### Getter
setMethod(
    "outputLibrary",
    "drugfindRCoreData",
    function(object) object@outputLibrary
)

#### Setter

setMethod(
    "setOutputLibrary",
    "drugfindRCoreData",
    function(object, value) {
        outputLibrary(object) <- value
    }
)

#### Replacer

setMethod(
    "outputLibrary<-",
    "drugfindRCoreData",
    function(object, value) {
        object@outputLibrary <- value
        validObject(object)
        object
    }
)

### -- Filtered Signature --

#### Generics
setGeneric(
    "filteredSignature",
    function(object) standardGeneric("filteredSignature")
)
setMethod(
    "filteredSignature",
    "drugfindRCoreData",
    function(object) object@filteredSignature
)

### -- Filter Thresholds --

#### Generics
setGeneric(
    "filterThresholdUp",
    function(object) standardGeneric("filterThresholdUp")
)

setGeneric(
    "setFilterThresholdUp",
    function(object, value) standardGeneric("setFilterThresholdUp")
)

setGeneric(
    "filterThresholdUp<-",
    function(object, value) standardGeneric("filterThresholdUp<-")
)

#### Getter
setMethod(
    "filterThresholdUp",
    "drugfindRCoreData",
    function(object) object@filterThresholdUp
)

#### Setter

setMethod(
    "setFilterThresholdUp",
    "drugfindRCoreData",
    function(object, value) {
        filterThresholdUp(object) <- value
    }
)

#### Replacer

setReplaceMethod(
    "filterThresholdUp",
    "drugfindRCoreData",
    function(object, value) {
        object@filterThresholdUp <- value
        validObject(object)
        object
    }
)


#### Generics

setGeneric(
    "filterThresholdDown",
    function(object) standardGeneric("filterThresholdDown")
)

setGeneric(
    "setFilterThresholdDown",
    function(object, value) standardGeneric("setFilterThresholdDown")
)

setGeneric(
    "filterThresholdDown<-",
    function(object, value) standardGeneric("filterThresholdDown<-")
)

#### Getter

setMethod(
    "filterThresholdDown",
    "drugfindRCoreData",
    function(object) object@filterThresholdDown
)

#### Setter

setMethod(
    "setFilterThresholdDown",
    "drugfindRCoreData",
    function(object, value) {
        filterThresholdDown(object) <- value
    }
)

#### Replacer

setReplaceMethod(
    "filterThresholdDown",
    "drugfindRCoreData",
    function(object, value) {
        object@filterThresholdDown <- value
        validObject(object)
        object
    }
)

### Generics

setGeneric(
    "filterThresholds",
    function(object) standardGeneric("filterThresholds")
)

setGeneric(
    "setFilterThresholds",
    function(object, value) standardGeneric("setFilterThresholds")
)

setGeneric(
    "filterThresholds<-",
    function(object, value) standardGeneric("filterThresholds<-")
)

#### Getter

setMethod(
    "filterThresholds",
    "drugfindRCoreData",
    function(object) {
        c(filterThresholdDown(object), filterThresholdUp(object))
    }
)

#### Setter

setMethod(
    "setFilterThresholds", "drugfindRCoreData",
    function(object, value) {
        filterThresholds(object) <- value
    }
)

#### Replacer

setReplaceMethod(
    "filterThresholds",
    "drugfindRCoreData",
    function(object, value) {
        if (length(value) == 1L) {
            object@filterThresholdDown <- -value
            object@filterThresholdUp <- value
        } else if (length(value) == 2L) {
            object@filterThresholdDown <- value[1L]
            object@filterThresholdUp <- value[2L]
        } else {
            stop("Incorrect number of items in length", call. = FALSE)
        }
        validObject(object)
        object
    }
)

### -- Concordance Limits --

#### Generics

setGeneric(
    "concordanceLimitUp",
    function(object) standardGeneric("concordanceLimitUp")
)
setGeneric(
    "setConcordanceLimitUp",
    function(object, value) standardGeneric("setConcordanceLimitUp")
)
setGeneric(
    "concordanceLimitUp<-",
    function(object, value) standardGeneric("concordanceLimitUp<-")
)

#### Getter

setMethod(
    "concordanceLimitUp",
    "drugfindRCoreData",
    function(object) object@concordanceLimitUp
)

#### Setter

setMethod(
    "setConcordanceLimitUp", "drugfindRCoreData",
    function(object, value) {
        concordanceLimitUp(object) <- value
    }
)

#### Replacer

setReplaceMethod(
    "concordanceLimitUp",
    "drugfindRCoreData",
    function(object, value) {
        object@concordanceLimitUp <- value
        validObject(object)
        object
    }
)

#### Generics

setGeneric(
    "concordanceLimitDown",
    function(object) standardGeneric("concordanceLimitDown")
)
setGeneric(
    "setConcordanceLimitDown",
    function(object, value) standardGeneric("setConcordanceLimitDown")
)
setGeneric(
    "concordanceLimitDown<-",
    function(object, value) standardGeneric("concordanceLimitDown<-")
)

#### Getter

setMethod(
    "concordanceLimitDown",
    "drugfindRCoreData",
    function(object) object@concordanceLimitDown
)

#### Setter

setMethod(
    "setConcordanceLimitDown",
    "drugfindRCoreData",
    function(object, value) {
        concordanceLimitDown(object) <- value
    }
)

#### Replacer

setReplaceMethod(
    "concordanceLimitDown",
    "drugfindRCoreData",
    function(object, value) {
        object@concordanceLimitDown <- value
        validObject(object)
        object
    }
)

#### Generics

setGeneric(
    "concordanceLimits",
    function(object) standardGeneric("concordanceLimits")
)

setGeneric(
    "setConcordanceLimits",
    function(object, value) standardGeneric("setConcordanceLimits")
)

setGeneric(
    "concordanceLimits<-",
    function(object, value) standardGeneric("concordanceLimits<-")
)

#### Getter

setMethod(
    "concordanceLimits",
    "drugfindRCoreData",
    function(object) {
        c(concordanceLimitDown(object), concordanceLimitUp(object))
    }
)

#### Setter

setMethod(
    "setConcordanceLimits",
    "drugfindRCoreData",
    function(object, value) {
        concordanceLimits(object) <- value
    }
)

#### Replacer

setReplaceMethod(
    "concordanceLimits",
    "drugfindRCoreData",
    function(object, value) {
        if (length(value) == 1L) {
            object@concordanceLimitDown <- -value
            object@concordanceLimitUp <- value
        } else if (length(value) == 2L) {
            object@concordanceLimitDown <- value[1L]
            object@concordanceLimitUp <- value[2L]
        } else {
            stop("Incorrect number of items in length", call. = FALSE)
        }
        validObject(object)
        object
    }
)

### -- Concordants --
setGeneric(
    "unfilteredConcordants",
    function(object) standardGeneric("unfilteredConcordants")
)
setMethod(
    "unfilteredConcordants",
    "drugfindRCoreData",
    function(object) object@unfilteredConcordants
)

setGeneric(
    "filteredConcordants",
    function(object) standardGeneric("filteredConcordants")
)
setMethod(
    "filteredConcordants",
    "drugfindRCoreData",
    function(object) object@filteredConcordants
)
