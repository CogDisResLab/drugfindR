#' @include utilities.R
#' @include getSignature.R prepareSignature.R
#' @include getConcordants.R consensusConcordants.R filterSignature.R
NULL

#' Investigate a Given Gene or Drug
#' `r lifecycle::badge("experimental")`
#'
#' This function takes the name of a gene or a drug and a
#' database to use to pull signatures
#' from and then queries iLINCS to get concordant signatures
#'
#' @param target The name of the gene or drug
#' @param inputLib One of "OE", "KD" or "CP". Marks the database to use.
#' @param outputLib One of "OE", "KD" or "CP". Marks the database to query.
#' @param filterThreshold The Filtering threshold.
#' @param similarityThreshold The Similarity Threshold
#' @param paired Logical. Whether to query iLINCS separately
#' for up and down regulated genes
#' @param inputCellLines A character vector of cell lines to
#' restrict our search for input signatures to.
#' @param outputCellLines A character vetor of cell lines to
#' restrict the output search to.
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
#' \donttest{
#'
#' # Search the whole iLINCS database for top concordant signatures for an
#' # ABL2 knockdown signature
#'
#' investigatedSignature <- investigateTarget("ABL2",
#'     inputLib = "KD",
#'     outputLib = "CP",
#'     filterThreshold = 0.5
#' )
#' }
#'
investigateTarget <- function(target,
                              inputLib, outputLib,
                              filterThreshold = 0.85,
                              similarityThreshold = 0.321,
                              paired = TRUE, inputCellLines = NULL,
                              outputCellLines = NULL, discordant = FALSE) {
    libs <- c("OE", "KD", "CP")

    if (!inputLib %in% libs || !outputLib %in% libs) {
        stop("Both input and output libraries must be one of 'OE', 'KD', 'CP'")
    }


    if (inputLib == "OE") {
        inputMetadata <- oeMetadata # nolint: object_usage_linter.
    } else if (inputLib == "KD") {
        inputMetadata <- kdMetadata # nolint: object_usage_linter.
    } else if (inputLib == "CP") {
        inputMetadata <- cpMetadata # nolint: object_usage_linter.
    } else {
        stop("Invalid inputLib")
    }


    if (!is.null(inputCellLines)) {
        filteredSignatureIds <- inputMetadata %>%
            dplyr::filter(
                stringr::str_to_lower(target) ==
                    stringr::str_to_lower(.data[["Source"]])
            ) %>%
            dplyr::filter(.data[["SourceCellLine"]] %in% inputCellLines) %>%
            dplyr::pull(.data[["SourceSignature"]])
    } else {
        filteredSignatureIds <- inputMetadata %>%
            dplyr::filter(
                stringr::str_to_lower(target) ==
                    stringr::str_to_lower(.data[["Source"]])
            ) %>%
            dplyr::pull(.data[["SourceSignature"]])
    }

    if (length(filteredSignatureIds) == 0L) {
        stop("No signatures match the given input criteria.")
    }

    allSignatures <- filteredSignatureIds %>%
        purrr::map(~ getSignature(.x))

    if (paired) {
        filteredUp <- allSignatures %>%
            purrr::map(~ filterSignature(.x,
                direction = "up",
                threshold = filterThreshold
            ))

        filteredDown <- allSignatures %>%
            purrr::map(~ filterSignature(.x,
                direction = "down",
                threshold = filterThreshold
            ))

        concordantUp <- filteredUp %>%
            purrr::map(~ getConcordants(.x, library = outputLib))

        concordantDown <- filteredDown %>%
            purrr::map(~ getConcordants(.x, library = outputLib))

        consensusTargets <- purrr::map2(
            concordantUp, concordantDown,
            ~ consensusConcordants(.x, .y,
                paired = paired,
                cellLine = outputCellLines,
                discordant = discordant,
                cutoff = similarityThreshold
            )
        )
    } else {
        filtered <- allSignatures %>%
            purrr::map(~ filterSignature(.x,
                direction = "any",
                threshold = filterThreshold
            ))

        concordants <- filtered %>%
            purrr::map(~ getConcordants(.x, library = outputLib))

        consensusTargets <- purrr::map(
            concordants,
            ~ consensusConcordants(.x,
                paired = paired,
                cellLine = outputCellLines,
                discordant = discordant,
                cutoff = similarityThreshold
            )
        )
    }

    augmented <- consensusTargets %>%
        purrr::map2(
            filteredSignatureIds,
            ~ dplyr::mutate(.x, SourceSignature = .y)
        ) %>%
        purrr::map_dfr(~ dplyr::inner_join(.x,
            inputMetadata,
            by = "SourceSignature"
        )) %>%
        dplyr::select(
            dplyr::any_of(c(
                "Source",
                "Target",
                "Similarity",
                "SourceSignature",
                "SourceCellLine",
                "SourceConcentration",
                "SourceTime",
                "TargetSignature",
                "TargetCellLine",
                "TargetConcentration",
                "TargetTime"
            ))
        )

    augmented
}
