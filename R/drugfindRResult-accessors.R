#' @include class-drugfindRResult.R
NULL


#' Accessor Methods for drugfindRResult Slots
#'
#' These functions provide internal access to the individual slots of a
#' `drugfindRResult` object used in the drugfindR pipeline.
#'
#' @param object A [`drugfindRResult`] object.
#' @param value The value to assign to a slot (for setters only).
#'
#' @return For getters, the value stored in the relevant slot.
#' For setters, the updated object.
#'
#' ## Getter methods
#'
#' - `coreDataObject(object)` returns the core gene signature.
#' - `sourceClass(object)` returns both thresholds as a vector.
#'
#' ## Setter methods
#'
#' - `coreDataObject(object)` returns the core gene signature.
#' - `sourceClass(object)` returns both thresholds as a vector.
#'
#' @name drugfindRResult-accessors
#' @aliases resultData sourceClass
#' resultData<- sourceClass<-
#'
#' @importFrom methods validObject
#'
#' @keywords internal

### -- Core Data Object --

#### Generics
setGeneric(
    "resultData",
    function(object, ...) standardGeneric("resultData")
)

setGeneric(
    "resultData<-",
    function(object, ...) standardGeneric("resultData<-")
)

#### Getter
setMethod(
    "resultData", "drugfindRResult",
    function(object, ...) object@result
)

#### Setter
setReplaceMethod(
    "resultData",
    "drugfindRResult",
    function(object, value) {
        object@result <- value
        validObject(object)
        object
    }
)

### -- Source Class --

#### Generics
setGeneric(
    "sourceClass",
    function(object, ...) standardGeneric("sourceClass")
)

setGeneric(
    "sourceClass<-",
    function(object, ...) standardGeneric("sourceClass<-")
)

#### Getter
setMethod(
    "sourceClass", "drugfindRResult",
    function(object, ...) object@sourceClass
)

#### Setter
setReplaceMethod(
    "sourceClass",
    "drugfindRResult",
    function(object, value) {
        object@sourceClass <- value
        validObject(object)
        object
    }
)
