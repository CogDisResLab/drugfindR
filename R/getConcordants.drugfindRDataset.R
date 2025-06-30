#' @include specialGenerics.R
#' @include drugfindRDataset-accessors.R
#' @include drugfindRResult.R
NULL

#' Get the concordant signatures from the designated library
#'
#' @param object An object of type [`drugfindRCoreData`]
#'
#' @return A copy of the input object, with the @unfilteredConcordants
#' slot populated with the concordants
#'
#' @importFrom httr2 request req_url_path_append req_url_query
#' @importFrom httr2 req_body_multipart req_method req_user_agent
#' @importFrom httr2 req_perform resp_status resp_body_string
#' @importFrom readr write_tsv
#' @importFrom purrr map flatten_dfr
#' @importFrom dplyr select mutate
#' @importFrom httr2 req_perform
#' @importFrom httr2 req_perform
#' @importFrom methods validObject
#'
#' @export
setMethod("getConcordants", "drugfindRDataset", function(object) {
    concardantsProcessed <- getConcordants(coreData(object))
    sourceClass <- sourceClass(object)
    result <- drugfindRResult(concardantsProcessed, sourceClass)
    result
})
