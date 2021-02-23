#' Filter the L1000 Signature
#'
#' This function filters the L1000 Signature to a given threshold, identifying
#' up-regulated or down-regulated or both up- and down-regulated genes
#'
#' @param signature A dataframe with the L1000 signature
#' @param direction Direction to filter to. Must be one of "up", "down" or "any". Defaults to "any"
#' @param threshold A Log Fold-Change Threshold to filter at. Is inclusive and defaults to 0.85
#'
#' @return a tibble with the filtered L1000 Signature
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
filter_signature <- function(signature, direction = "any", threshold = 0.85) {
  stopifnot("data.frame" %in% class(signature))

  if (!direction %in% c("up", "down", "any")) {
    stop("Direction must be one of 'up', 'down' or 'any'")
  }

  if (direction == "up") {
    filtered <- signature %>%
      dplyr::filter(.data$Value_LogDiffExp >= threshold)
  } else if (direction == "down") {
    filtered <- signature %>%
      dplyr::filter(.data$Value_LogDiffExp <= -threshold)
  } else {
    filtered <- signature %>%
      dplyr::filter(abs(.data$Value_LogDiffExp) >= threshold)
  }

  filtered
}
