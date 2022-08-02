#' Get Concordant Signatures from iLINCS
#'
#' This function takes a full or filtered signature and gets concordant signatures
#' from any of the 3 LINCS databases in iLINCS. This can get Overexpression,
#' Knockdown or Chemical Perturbagen signatures.
#'
#' @param signature A data frame with the names of genes, their expression value
#' and optionally their p-value
#' @param library The Library you want to search. Must be one of "OE", "KD" or "CP"
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
#'
#' @examples
#' TRUE
get_concordants <- function(signature, library = "CP") {
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

  if (!library %in% c("OE", "KD", "CP")) {
    stop("library must be one of 'OE', 'KD' or 'CP'")
  }

  url <- "http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze"
  query <- list(lib = lib_map[library])
  body <- list(file = httr::upload_file(signature_file))

  request <- httr::POST(url, query = query, body = body)

  if (httr::status_code(request) == 200) {
    concordants <- httr::content(request) %>%
      purrr::map("concordanceTable") %>%
      purrr::flatten_dfr() %>%
      dplyr::select(dplyr::any_of(c(
        "signatureid", "compound", "treatment",
        "concentration", "time", "cellline", "similarity", "pValue"
      ))) %>%
      dplyr::mutate(
        similarity = round(.data$similarity, 8),
        pValue = round(.data$pValue, 20)
      )
  } else if (library %in% c("OE", "KD")) {
    concordants <- tibble::tibble(
      signatureid = NA,
      treatment = NA,
      time = NA,
      cellline = NA,
      similarity = NA,
      pValue = NA
    ) %>%
      dplyr::filter(!is.na(.data$signatureid))
  } else {
    concordants <- tibble::tibble(
      signatureid = NA,
      compound = NA,
      concentration = NA,
      time = NA,
      cellline = NA,
      similarity = NA,
      pValue = NA
    ) %>%
      dplyr::filter(!is.na(.data$signatureid))
  }

  concordants
}
