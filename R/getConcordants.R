#' Get Concordant Signatures from iLINCS
#' `r lifecycle::badge("experimental")`
#'
#' This function takes a full or filtered signature
#' and gets concordant signatures
#' from any of the 3 LINCS databases in iLINCS. This can get Overexpression,
#' Knockdown or Chemical Perturbagen signatures.
#'
#' @param signature A data frame with the names of genes, their expression value
#' and optionally their p-value
#' @param ilincsLibrary The Library you want to search.
#' Must be one of "OE", "KD" or "CP"
#' for Overexpression, Knockdown or Chemical Perturbagens
#'
#' @return A tibble with the list of concordant and discordant signatures
#' @export
#'
#' @importFrom readr write_tsv
#' @importFrom httr POST status_code content upload_file
#' @importFrom httr2 request req_method req_url_query req_user_agent req_url_path_append req_perform resp_status resp_body_json resp_body_string
#' @importFrom purrr map flatten_dfr
#' @importFrom dplyr select any_of mutate filter
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stringr str_glue
#' @importFrom curl form_file
#'
#' @examples
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
#'
#' # Get concordant gene knockdown signatures
#'
#' concordant_signatures <- getConcordants(kdSignature, ilincsLibrary = "KD")
#'
#' head(concordant_signatures)
getConcordants <- function(signature, ilincsLibrary = "CP") {
    if (!"data.frame" %in% class(signature)) {
        stop("signature must be a data frame or data frame like object")
    } else {
        signatureFile <- tempfile(pattern = "ilincs_sig", fileext = ".xls")
        signature |>
            readr::write_tsv(signatureFile)
    }

    sigDirection <- if (all(signature[["Value_LogDiffExp"]] > 0L)) {
        "Up"
    } else if (all(signature[["Value_LogDiffExp"]] < 0L)) {
        "Down"
    } else {
        "Any"
    }

    libMap <- c(
        OE = "LIB_11",
        KD = "LIB_6",
        CP = "LIB_5"
    )

    stopIfInvalidLibraries(ilincsLibrary)

    request <- httr2::request(.ilincsBaseUrl()) |>
        httr2::req_url_path_append("SignatureMeta") |>
        httr2::req_url_path_append("uploadAndAnalyze") |>
        httr2::req_url_query(lib = libMap[ilincsLibrary]) |>
        httr2::req_body_multipart(file = curl::form_file(signatureFile)) |>
        httr2::req_method("POST") |>
        httr2::req_user_agent(stringr::str_glue("drugfindR/v{packageVersion('drugfindR')}; https://github.com/CogDisResLab/drugfindR")) |>
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
        return(concordants)
    } else {
        stop(httr2::resp_status(request), " ", httr2::resp_body_string(request))
    }
}
