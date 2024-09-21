#' Get the L1000 Signature from iLINCS
#' `r lifecycle::badge("stable")`
#'
#' This function acts as the entrypoint to the iLINCS database.
#' This takes in an ID and returns the signature after making a
#' call to the iLINCS
#' database. The default mode for `drugfindR` is to use L1000
#' signatures. However,
#' if you are trying to retrieve a different transcriptomic signature,
#' that is also supported
#' by setting the `l1000` parameter to `FALSE`.
#'
#' @param sigId character. The ilincs signature_id
#' @param l1000 boolean. If you have a known l1000 signature
#'
#' @return a tibble with the L1000 Signature
#' @export
#'
#' @importFrom httr2 request req_method req_url_query req_user_agent req_url_path_append req_perform resp_status resp_body_json resp_body_string
#' @importFrom tibble tibble as_tibble
#' @importFrom rlang .data
#' @importFrom dplyr select
#' @importFrom purrr map_dfr
#' @importFrom S4Vectors DataFrame
#'
#' @examples
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
getSignature <- function(sigId, l1000 = TRUE) {
    if (l1000) {
        numGenes <- 978L
    } else {
        numGenes <- Inf
    }

    request <- httr2::request(.ilincsBaseUrl()) |>
        httr2::req_url_path_append("ilincsR") |>
        httr2::req_url_path_append("downloadSignature") |>
        httr2::req_url_query(sigID = sigId, noOfTopGenes = numGenes) |>
        httr2::req_method("POST") |>
        httr2::req_user_agent("drugfindR/v0.99.980; https://github.com/CogDisResLab/drugfindR") |>
        httr2::req_perform()


    if (httr2::resp_status(request) == 200L) {
        signature <- httr2::resp_body_json(request) |>
            purrr::map("signature") |>
            purrr::flatten_dfr() |>
            dplyr::select(-"PROBE") |>
            dplyr::mutate(
                ID_geneid = as.character(.data[["ID_geneid"]]),
                Value_LogDiffExp = round(.data[["Value_LogDiffExp"]], 12L),
                Significance_pvalue = round(.data[["Significance_pvalue"]], 12L)
            )
    } else {
        stop(httr2::resp_status(request), " ", httr2::resp_body_string(request))
    }
}
