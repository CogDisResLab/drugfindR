#' Get the L1000 Signature from iLINCS
#'
#' `r lifecycle::badge("experimental")`
#'
#' This function acts as the entrypoint to the iLINCS database.
#' This takes in an ID and returns the signature after making a call to the iLINCS
#' database. The default mode for `drugfindR` is to use L1000 signatures. However,
#' if you are trying to retrieve a different transcriptomic signature, that is also supported
#' by setting the `l1000` parameter to `FALSE`.
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
#' \dontrun{
#'
#' # Get the L1000 signature for LINCSKD_28
#' kd_signature <- getSignature("LINCSKD_28")
#'
#' # Get the non-L1000 signature for EBI_1001
#' ebi_signature <- getSignature("EBI_1001")
#' }
getSignature <- function(sigId, l1000 = TRUE) {
    url <- "http://www.ilincs.org/api/ilincsR/downloadSignature"

    if (l1000) {
        numGenes <- 978L
    } else {
        numGenes <- 25000L
    }

    query <- list(sigID = sigId, noOfTopGenes = numGenes)

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
        stop(httr::status_code(request), " ", httr::content(request))
    }
}
