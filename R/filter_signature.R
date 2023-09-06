#' Filter the L1000 Signature
#'
#' `r lifecycle::badge("experimental")`
#'
#' This function filters the L1000 Signature to a given threshold, identifying
#' up-regulated or down-regulated or both up- and down-regulated genes
#'
#' @param signature A dataframe with the L1000 signature
#' @param direction Direction to filter to. Must be one of "up", "down" or "any". Defaults to "any"
#' @param threshold A Log Fold-Change Threshold to filter at. Cannot be specified with prop
#' @param prop A proportion of genes to take from top and bottom. Cannot be specified with threshold
#'
#' @return a tibble with the filtered L1000 Signature
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom stats quantile
#'
#' @examples
#' TRUE
filter_signature <- function(signature, direction = "any", threshold = NULL, prop = NULL) {
    stopifnot("data.frame" %in% class(signature))

    if (!is.null(threshold) && !is.null(prop)) {
        stop("Only one of prop or threshold can be specified")
    } else if (is.null(threshold) && is.null(prop)) {
        stop("One of prop or threshold must be specified")
    }

    if (!direction %in% c("up", "down", "any")) {
        stop("Direction must be one of 'up', 'down' or 'any'")
    }

    if (!is.null(threshold)) {
        if (length(threshold) == 2L) {
            down_threshold <- threshold[[1L]]
            up_threshold <- threshold[[2L]]
        } else if (length(threshold) == 1L) {
            down_threshold <- -threshold
            up_threshold <- threshold
        } else {
            stop("Threshold must be specified as one or two values")
        }
    } else if (!is.null(prop)) {
        down_threshold <- quantile(signature[["Value_LogDiffExp"]], prop)
        up_threshold <- quantile(signature[["Value_LogDiffExp"]], 1.0 - prop)
    }

    if (direction == "up") {
        filtered <- signature %>%
            dplyr::filter(.data[["Value_LogDiffExp"]] >= up_threshold)
    } else if (direction == "down") {
        filtered <- signature %>%
            dplyr::filter(.data[["Value_LogDiffExp"]] <= down_threshold)
    } else {
        filtered <- signature %>%
            dplyr::filter(
                .data[["Value_LogDiffExp"]] >= up_threshold |
                    .data[["Value_LogDiffExp"]] <= down_threshold
            )
    }

    filtered
}
