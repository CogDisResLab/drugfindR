#' Get the L1000 Signature from iLINCS
#'
#' @param sig_id character. The ilincs signature_id
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
get_signature <- function(sig_id) {
  url <- "http://www.ilincs.org/api/ilincsR/downloadSignature"
  query = list(sigID = sig_id, noOfTopGenes = 978)

  request <- httr::POST(url, query = query)

  if (httr::status_code(request) == 200) {
    signature <- httr::content(request) %>%
      purrr::map("signature") %>%
      purrr::flatten_dfr() %>%
      dplyr::select(-.data$PROBE) %>%
      dplyr::mutate(Value_LogDiffExp = round(Value_LogDiffExp, 12),
                    Significance_pvalue = round(Significance_pvalue, 12))
  } else {
    signature <- tibble::tibble(
      signatureID = rep(NA, 978),
      ID_geneid = rep(NA, 978),
      Name_GeneSymbol = rep(NA, 978),
      Value_LogDiffExp = rep(NA, 978),
      Significance_pvalue = rep(NA, 978)
    )
  }

  signature
}
