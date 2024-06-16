#' Filter the L1000 Signature
#' `r lifecycle::badge("experimental")`
#'
#' This function filters the L1000 Signature to a given threshold, identifying
#' up-regulated or down-regulated or both up- and down-regulated genes
#'
#' @param signature A dataframe with the L1000 signature
#' @param direction Direction to filter to. Must be one of "up",
#' "down" or "any". Defaults to "any"
#' @param threshold A Log Fold-Change Threshold to filter at.
#' This can either be a single value or a vector of two values.
#' If a single value is given, then it is assumed to be a symmetric threshold.
#' If two values are given, then the first value is the down-regulated
#' threshold and the second value is
#' the up-regulated threshold. Cannot be specified with prop
#' @param prop A proportion of genes to take from top and bottom.
#' Cannot be specified with threshold
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
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
#'
#' # Filter signature by a specific threshold
#' filteredSignature <- filterSignature(kdSignature, threshold = 0.5)
#'
#' # Filter signature by a proportion
#'
#' filteredSignature <- filterSignature(kdSignature, prop = 0.1)
#'
#' # Filter Signature to up-regulated genes only by a threshold
#'
#' filteredSignature <- filterSignature(kdSignature,
#'     direction = "up", threshold = 0.5
#' )
#'
#' # Filter the signature using differing thresholds for up and
#' # down-regulated genes
#'
#' filteredSignature <- filterSignature(kdSignature,
#'     threshold = c(-0.75, 0.5)
#' )
filterSignature <- function(
    signature, direction = "any",
    threshold = NULL, prop = NULL) {
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
            downThreshold <- threshold[[1L]]
            upThreshold <- threshold[[2L]]
        } else if (length(threshold) == 1L) {
            downThreshold <- -threshold
            upThreshold <- threshold
        } else {
            stop("Threshold must be specified as one or two values")
        }
    } else if (!is.null(prop)) {
        downThreshold <- quantile(signature[["Value_LogDiffExp"]], prop)
        upThreshold <- quantile(signature[["Value_LogDiffExp"]], 1.0 - prop)
    }

    if (direction == "up") {
        filtered <- signature %>%
            dplyr::filter(.data[["Value_LogDiffExp"]] >= upThreshold)
    } else if (direction == "down") {
        filtered <- signature %>%
            dplyr::filter(.data[["Value_LogDiffExp"]] <= downThreshold)
    } else {
        filtered <- signature %>%
            dplyr::filter(
                .data[["Value_LogDiffExp"]] >= upThreshold |
                    .data[["Value_LogDiffExp"]] <= downThreshold
            )
    }

    filtered
}
