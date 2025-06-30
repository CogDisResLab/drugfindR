#' @include utilities.R
NULL


### Method to populate filtered signature slot
setGeneric(
    "filterSignature",
    function(object) standardGeneric("filterSignature")
)


### Method for populating the unfiltered concordants slot
setGeneric(
    "getConcordants",
    function(object) standardGeneric("getConcordants")
)

### Method to populate filtered concordants slot
setGeneric(
    "consensusConcordants",
    function(object) standardGeneric("consensusConcordants")
)
