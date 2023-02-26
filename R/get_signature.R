#' Get the L1000 Signature from iLINCS
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
        num_genes <- 978
    } else {
        num_genes <- 25000
    }

    query <- list(sigID = sig_id, noOfTopGenes = num_genes)

    request <- httr::POST(url, query = query)

    if (httr::status_code(request) == 200) {
        signature <- httr::content(request) %>%
            purrr::map("signature") %>%
            purrr::flatten_dfr() %>%
            dplyr::select(-PROBE) %>%
            dplyr::mutate(
                Value_LogDiffExp = round(.data$Value_LogDiffExp, 12),
                Significance_pvalue = round(.data$Significance_pvalue, 12)
            )
    } else {
        signature <- tibble::tibble(
            signatureID = rep(NA, num_genes),
            ID_geneid = rep(NA, num_genes),
            Name_GeneSymbol = rep(NA, num_genes),
            Value_LogDiffExp = rep(NA, num_genes),
            Significance_pvalue = rep(NA, num_genes)
        )
    }

    signature
}
