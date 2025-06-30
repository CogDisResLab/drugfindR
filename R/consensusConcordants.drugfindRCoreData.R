#' @include specialGenerics.R
#' @include class-drugfindRCoreData.R
NULL

#' Generate a Consensus list of Targets
#' `r lifecycle::badge("stable")`
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
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
#'
#' # Get concordant gene knockdown signatures
#' concordantSignatures <- getConcordants(kdSignature, ilincsLibrary = "KD")
#'
#' # Get the consensus list of signatures with defaults
#' consensus <- consensusConcordants(concordantSignatures)
#'
#' # Get the consensus list of signatures with a different cutoff
#' consensus <- consensusConcordants(concordantSignatures,
#'     cutoff = 0.5
#' )
#'
#' # Get the consensus list of signatures with a specified cell lines
#' consensus <- consensusConcordants(concordantSignatures,
#'     cellLine = c("A549", "MCF7")
#' )
#'
#' # Doing a paired analysis
#' filteredUp <- filterSignature(kdSignature,
#'     direction = "up", threshold = 0.5
#' )
#' filteredDown <- filterSignature(kdSignature,
#'     direction = "down", threshold = -0.5
#' )
#'
#' concordants_up <- getConcordants(filteredUp, ilincsLibrary = "KD")
#' concordants_down <- getConcordants(filteredDown, ilincsLibrary = "KD")
#'
#' consensus <- consensusConcordants(concordants_up,
#'     concordants_down,
#'     paired = TRUE
#' )
#'
setMethod(
    "consensusConcordants", "drugfindRCoreData",
    function(object) {
        if (!is.null(object@cellLines)) {
            concordants <- dplyr::filter(
                object@unfilteredConcordants,
                cellline %in% object@cellLines
            )
        }

        filtered <- concordants |>
            dplyr::filter(
                similarity < concordanceLimitDown(object) |
                    similarity > concordanceLimitUp(object)
            ) |>
            dplyr::group_by(
                dplyr::across(dplyr::any_of(c("treatment", "compound")))
            ) |>
            dplyr::filter(
                abs(similarity) == max(abs(similarity))
            ) |>
            dplyr::select(
                dplyr::any_of(c(
                    "signatureid", "treatment", "compound", "cellline", "time",
                    "concentration", "similarity", "sig_direction", "pValue"
                ))
            ) |>
            dplyr::arrange(dplyr::desc(abs(similarity))) |>
            dplyr::rename_with(targetRename) |>
            dplyr::ungroup()

        object@filteredConcordants <- filtered
        validObject(object)
        object
    }
)
