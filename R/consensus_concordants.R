#' Rename the Target-Related Columns
#'
#' @param input_names A character vector of input_names
#'
#' @return A character vector of new names
#'
#' @examples
#' TRUE
target_rename <- function(input_names) {
  if ("treatment" %in% input_names) {
    new_cols <- c("TargetSignature", "Target", "TargetCellLine",
                  "TargetTime", "Similarity")
  } else {
    new_cols <- c("TargetSignature", "Target", "TargetCellLine",
                  "TargetTime", "TargetConcentration", "Similarity")
  }

  new_cols
}

#' Generate a Consensus list of Targets
#'
#' This function takes a list of (optionally split) concordance dataframes and returns
#' a ranked list of gene or drug targets that have been chose for their maximal
#' similarity to the signature
#'
#' @param ... One or Two (see paired) Data Frames with the concordants
#' @param paired Logical indicating whether you split the dataframes by up and down regulated in prior analysis
#' @param cutoff A similarity cutoff value. Defaults to 0.321
#' @param cell_line A character vector of Cell Lines you are interested in.
#' @param discordant Logical indicating whether you want discordant signatures
#'
#' @return A tibble with the filtered and deduplicated results
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange any_of group_by across select bind_rows rename_with ungroup
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
consensus_concordants <- function(..., paired = FALSE, cutoff = 0.321,
                                  cell_line = NULL, discordant = FALSE) {

  if (paired & length(list(...)) != 2) {
    stop("Paired analysis requires two data frames")
  } else if (!paired & length(list(...)) != 1) {
    stop("Unpaired analysis requires only one dataframe")
  }

  concordants <- dplyr::bind_rows(list(...))

  if (!is.null(cell_line)) {
    concordants <- concordants %>%
      dplyr::filter(.data$cellline %in% cell_line)
  }

  if (discordant) {
    filtered <- concordants %>%
      dplyr::filter(.data$similarity <= -cutoff) %>%
      dplyr::group_by(
        dplyr::across(dplyr::any_of(c("treatment", "compound")))
      ) %>%
      dplyr::filter(abs(.data$similarity) == max(abs(.data$similarity))) %>%
      dplyr::select(.data$signatureid, dplyr::any_of(c("treatment", "compound")),
                    .data$cellline, .data$time, dplyr::any_of(c("concentration")),
                    .data$similarity) %>%
      dplyr::arrange(dplyr::desc(abs(.data$similarity))) %>%
      dplyr::rename_with(target_rename) %>%
      dplyr::ungroup()
  } else {
    filtered <- concordants %>%
      dplyr::filter(.data$similarity >= cutoff) %>%
      dplyr::group_by(
        dplyr::across(dplyr::any_of(c("treatment", "compound")))
      ) %>%
      dplyr::filter(abs(.data$similarity) == max(abs(.data$similarity))) %>%
      dplyr::select(.data$signatureid, dplyr::any_of(c("treatment", "compound")),
                    .data$cellline, .data$time, dplyr::any_of(c("concentration")),
                    .data$similarity) %>%
      dplyr::arrange(dplyr::desc(abs(.data$similarity))) %>%
      dplyr::rename_with(target_rename) %>%
      dplyr::ungroup()
  }

  filtered
}
