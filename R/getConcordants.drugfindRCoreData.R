#' @include utilities.R
NULL

setGeneric(
    "getConcordants",
    function(object) standardGeneric("getConcordants")
)

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
setMethod("getConcordants", "drugfindRCoreData", function(object) {
    signatureFile <- tempfile(pattern = "ilincs_sig", fileext = ".xls")

    readr::write_tsv(object@signature, signatureFile)

    sigDirection <- determineSignatureDirection(object@signature)

    request <- httr2::request(.ilincsBaseUrl()) |>
        httr2::req_url_path_append("SignatureMeta") |>
        httr2::req_url_path_append("uploadAndAnalyze") |>
        httr2::req_url_query(lib = .return_library(object@outputLibrary)) |>
        httr2::req_body_multipart(file = curl::form_file(signatureFile)) |>
        httr2::req_method("POST") |>
        httr2::req_user_agent(.returnUserAgent()) |>
        httr2::req_perform()

    if (httr2::resp_status(request) == 200L) {
        concordants <- httr2::resp_body_json(request) |>
            purrr::map("concordanceTable") |>
            purrr::flatten_dfr() |>
            dplyr::select(dplyr::any_of(c(
                "signatureid", "compound", "treatment",
                "concentration", "time", "cellline", "similarity", "pValue"
            ))) |>
            dplyr::mutate(
                similarity = round(.data[["similarity"]], 8L),
                pValue = round(.data[["pValue"]], 20L),
                sig_direction = sigDirection
            )
        object@unfilteredConcordants <- concordants
        validObject(object)
        object
    } else {
        stop(
            httr2::resp_status(request), " ",
            httr2::resp_body_string(request),
            call. = FALSE
        )
    }
})
