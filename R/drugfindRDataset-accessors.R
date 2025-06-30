#' @include class-drugfindRDataset.R
NULL


#' Accessor Methods for drugfindRDataset Slots
#'
#' These functions provide internal access to the individual slots of a
#' `drugfindRDataset` object used in the drugfindR pipeline.
#'
#' @param object A [`drugfindRDataset`] object.
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
#' @name drugfindRDataset-accessors
#' @aliases coreData sourceClass
#' coreData<- sourceClass<-
#'
#' @importFrom methods validObject
#'
#' @keywords internal

### -- Core Data Object --

#### Generics
setGeneric(
    "coreData",
    function(object, ...) standardGeneric("coreData")
)

setGeneric(
    "coreData<-",
    function(object, ...) standardGeneric("coreData<-")
)

#### Getter
setMethod(
    "coreData", "drugfindRDataset",
    function(object, ...) object@core
)

#### Setter
setReplaceMethod(
    "coreData",
    "drugfindRDataset",
    function(object, value) {
        object@core <- value
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
    "sourceClass", "drugfindRDataset",
    function(object, ...) object@sourceClass
)

#### Setter
setReplaceMethod(
    "sourceClass",
    "drugfindRDataset",
    function(object, value) {
        object@sourceClass <- value
        validObject(object)
        object
    }
)
