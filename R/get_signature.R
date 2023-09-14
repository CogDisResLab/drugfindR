#' Get the L1000 Signature from iLINCS
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param sig_id character. The ilincs signature_id
#' @param l1000 boolean. If you have a known l1000 signature
#'
#' @return a tibble with the L1000 Signature
#' @export
#'
#' @importFrom httr POST content status_code
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr select
#' @importFrom purrr map_dfr
#'
#' @examples
#' TRUE
get_signature <- function(sig_id, l1000 = TRUE) {
    url <- "http://www.ilincs.org/api/ilincsR/downloadSignature"

    if (l1000) {
        num_genes <- 978L
    } else {
        num_genes <- 25000L
    }

    query <- list(sigID = sig_id, noOfTopGenes = num_genes)

    request <- httr::POST(url, query = query)

    if (httr::status_code(request) == 200L) {
        signature <- httr::content(request) %>%
            purrr::map("signature") %>%
            purrr::flatten_dfr() %>%
            dplyr::select(-"PROBE") %>%
            dplyr::mutate(
                ID_geneid = as.character(.data[["ID_geneid"]]),
                Value_LogDiffExp = round(.data[["Value_LogDiffExp"]], 12L),
                Significance_pvalue = round(.data[["Significance_pvalue"]], 12L)
            )
        signature
    } else {
        stop("Error: ", httr::status_code(request))
    }
}
