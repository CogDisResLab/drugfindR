#' Rename the Target-Related Columns
#'
#' @param inputNames A character vector of input_names
#'
#' @return A character vector of new names
#'
#' @examples
#' TRUE
targetRename <- function(inputNames) {
    if ("treatment" %in% inputNames) {
        newCols <- c(
            "TargetSignature", "Target", "TargetCellLine",
            "TargetTime", "Similarity", "SignatureDirection"
        )
    } else {
        newCols <- c(
            "TargetSignature", "Target", "TargetCellLine",
            "TargetTime", "TargetConcentration", "Similarity",
            "SignatureDirection"
        )
    }

    newCols
}

#' Generate a Consensus list of Targets
#'
#' `r lifecycle::badge("experimental")`
#'
#' This function takes a list of (optionally split)
#' concordance dataframes and returns
#' a ranked list of gene or drug targets that have
#' been chose for their maximal
#' similarity to the signature
#'
#' @param ... One or Two (see paired) Data Frames with the concordants
#' @param paired Logical indicating whether you split the
#' dataframes by up and down regulated in prior analysis
#' @param cutoff A similarity cutoff value. Defaults to 0.321
#' @param cellLine A character vector of Cell Lines you are interested in.
#'
#' @return A tibble with the filtered and deduplicated results
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange any_of group_by across
#' select bind_rows rename_with ungroup
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
consensusConcordants <- function(...,
                                 paired = FALSE,
                                 cutoff = 0.321,
                                 cellLine = NULL) {
    dots <- list(...)
    if (paired && length(dots) != 2L) {
        stop("Paired analysis requires two data frames")
    } else if (!paired && length(dots) != 1L) {
        stop("Unpaired analysis requires only one dataframe")
    }

    concordants <- dplyr::bind_rows(dots)

    if (!is.null(cellLine)) {
        concordants <- concordants %>%
            dplyr::filter(.data[["cellline"]] %in% cellLine)
    }

    filtered <- concordants %>%
        dplyr::filter(abs(.data[["similarity"]]) >= cutoff) %>%
        dplyr::group_by(
            dplyr::across(dplyr::any_of(c("treatment", "compound")))
        ) %>%
        dplyr::filter(
            abs(.data[["similarity"]]) == max(abs(.data[["similarity"]]))
        ) %>%
        dplyr::select(
            dplyr::any_of(c(
                "signatureid", "treatment", "compound", "cellline", "time",
                "concentration", "similarity", "sig_direction"
            ))
        ) %>%
        dplyr::arrange(dplyr::desc(abs(.data[["similarity"]]))) %>%
        dplyr::rename_with(targetRename) %>%
        dplyr::ungroup()

    filtered
}
