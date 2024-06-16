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
#' @importFrom purrr map flatten_dfr
#' @importFrom dplyr select any_of mutate filter
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom magrittr %>%
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
        signature %>%
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

    if (!ilincsLibrary %in% c("OE", "KD", "CP")) {
        stop("library must be one of 'OE', 'KD' or 'CP'")
    }

    url <- "http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze"
    query <- list(lib = libMap[ilincsLibrary])
    body <- list(file = httr::upload_file(signatureFile))

    request <- httr::POST(url, query = query, body = body)

    if (httr::status_code(request) == 200L) {
        concordants <- httr::content(request) %>%
            purrr::map("concordanceTable") %>%
            purrr::flatten_dfr() %>%
            dplyr::select(dplyr::any_of(c(
                "signatureid", "compound", "treatment",
                "concentration", "time", "cellline", "similarity", "pValue"
            ))) %>%
            dplyr::mutate(
                similarity = round(.data[["similarity"]], 8L),
                pValue = round(.data[["pValue"]], 20L),
                sig_direction = sigDirection
            )
        return(concordants)
    } else {
        httr::stop_for_status(request, "get concardant signatures")
    }
}
