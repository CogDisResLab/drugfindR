#' Get Concordant Signatures from iLINCS
#'
#' `r lifecycle::badge("experimental")`
#'
#' This function takes a full or filtered signature and gets concordant signatures
#' from any of the 3 LINCS databases in iLINCS. This can get Overexpression,
#' Knockdown or Chemical Perturbagen signatures.
#'
#' @param signature A data frame with the names of genes, their expression value
#' and optionally their p-value
#' @param ilincs_library The Library you want to search. Must be one of "OE", "KD" or "CP"
#' @param sig_direction The direction of the signature. Must be one of "Up" or "Down"
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
#' TRUE
get_concordants <- function(signature, ilincs_library = "CP", sig_direction = NULL) {
    if (!"data.frame" %in% class(signature)) {
        stop("signature must be a data frame or data frame like object")
    } else {
        signature_file <- tempfile(pattern = "ilincs_sig", fileext = ".xls")
        signature %>%
            readr::write_tsv(signature_file)
    }

    lib_map <- c(
        OE = "LIB_11",
        KD = "LIB_6",
        CP = "LIB_5"
    )

    if (!ilincs_library %in% c("OE", "KD", "CP")) {
        stop("library must be one of 'OE', 'KD' or 'CP'")
    }

    url <- "http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze"
    query <- list(lib = lib_map[ilincs_library])
    body <- list(file = httr::upload_file(signature_file))

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
                sig_direction = sig_direction
            )
        return(concordants)
    } else {
        httr::stop_for_status(request, "get concardant signatures")
    }
}
