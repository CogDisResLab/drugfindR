#' Investigate a Given Gene or Drug
#'
#' `r lifecycle::badge("experimental")`
#'
#' This function takes the name of a gene or a drug and a database to use to pull signatures
#' from and then queries iLINCS to get concordant signatures
#'
#' @param target The name of the gene or drug
#' @param input_lib One of "OE", "KD" or "CP". Marks the database to use.
#' @param output_lib One of "OE", "KD" or "CP". Marks the database to query.
#' @param filter_threshold The Filtering threshold.
#' @param similarity_threshold The Similarity Threshold
#' @param paired Logical. Whether to query iLINCS separately for up and down regulated genes
#' @param input_cell_lines A character vector of cell lines to restrict our search for input signatures to.
#' @param output_cell_lines A character vetor of cell lines to restrict the output search to.
#' @param discordant Logical. Whether to look for discordant signatures
#'
#' @return A tibble with the the similarity scores and signature metadata
#' @export
#'
#' @importFrom dplyr filter pull select any_of inner_join
#' @importFrom stringr str_to_lower
#' @importFrom purrr map map2 map_dfr
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
#' TRUE
investigate_target <- function(
        target, input_lib, output_lib,
        filter_threshold = 0.85, similarity_threshold = 0.321,
        paired = TRUE, input_cell_lines = NULL,
        output_cell_lines = NULL, discordant = FALSE) {
    libs <- c("OE", "KD", "CP")

    if (!input_lib %in% libs || !output_lib %in% libs) {
        stop("Both input and output libraries must be one of 'OE', 'KD', 'CP'")
    }

    if (missing(input_lib) || missing(output_lib)) {
        stop("Please specify both input and output libraries")
    }

    if (input_lib == "OE") {
        input_metadata <- oe_metadata # nolint: object_usage_linter.
    } else if (input_lib == "KD") {
        input_metadata <- kd_metadata # nolint: object_usage_linter.
    } else if (input_lib == "CP") {
        input_metadata <- cp_metadata # nolint: object_usage_linter.
    } else {
        stop("Invalid input_lib")
    }


    if (!is.null(input_cell_lines)) {
        filtered_signature_ids <- input_metadata %>%
            dplyr::filter(stringr::str_to_lower(target) == stringr::str_to_lower(.data[["Source"]])) %>%
            dplyr::filter(.data[["SourceCellLine"]] %in% input_cell_lines) %>%
            dplyr::pull(.data[["SourceSignature"]])
    } else {
        filtered_signature_ids <- input_metadata %>%
            dplyr::filter(stringr::str_to_lower(target) == stringr::str_to_lower(.data[["Source"]])) %>%
            dplyr::pull(.data[["SourceSignature"]])
    }

    if (length(filtered_signature_ids) == 0L) {
        stop("No signatures match the given input criteria.")
    }

    all_signatures <- filtered_signature_ids %>%
        purrr::map(~ get_signature(.x))

    if (paired) {
        filtered_up <- all_signatures %>%
            purrr::map(~ filter_signature(.x, direction = "up", threshold = filter_threshold))

        filtered_down <- all_signatures %>%
            purrr::map(~ filter_signature(.x, direction = "down", threshold = filter_threshold))

        concordant_up <- filtered_up %>%
            purrr::map(~ get_concordants(.x, library = output_lib))

        concordant_down <- filtered_down %>%
            purrr::map(~ get_concordants(.x, library = output_lib))

        consensus_targets <- purrr::map2(
            concordant_up, concordant_down,
            ~ consensus_concordants(.x, .y,
                paired = paired,
                cell_line = output_cell_lines,
                discordant = discordant,
                cutoff = similarity_threshold
            )
        )
    } else {
        filtered <- all_signatures %>%
            purrr::map(~ filter_signature(.x, direction = "any", threshold = filter_threshold))

        concordants <- filtered %>%
            purrr::map(~ get_concordants(.x, library = output_lib))

        consensus_targets <- purrr::map(
            concordants,
            ~ consensus_concordants(.x,
                paired = paired,
                cell_line = output_cell_lines,
                discordant = discordant,
                cutoff = similarity_threshold
            )
        )
    }

    augmented <- consensus_targets %>%
        purrr::map2(filtered_signature_ids, ~ dplyr::mutate(.x, SourceSignature = .y)) %>%
        purrr::map_dfr(~ dplyr::inner_join(.x, input_metadata, by = "SourceSignature")) %>%
        dplyr::select(
            dplyr::any_of(c(
                "Source", "Target", "Similarity", "SourceSignature",
                "SourceCellLine", "SourceConcentration", "SourceTime", "TargetSignature",
                "TargetCellLine", "TargetConcentration", "TargetTime"
            ))
        )

    augmented
}
