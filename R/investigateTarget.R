#' @include utilities.R
#' @include getSignature.R prepareSignature.R
#' @include getConcordants.R consensusConcordants.R filterSignature.R
NULL

#' Investigate concordant signatures for a gene or drug
#' `r lifecycle::badge("stable")`
#'
#' Given the name of a target (gene knockdown/overexpression or compound)
#' this high–level convenience wrapper:
#'
#' 1. Locates iLINCS source signatures for the target in the specified input library.
#' 1. Optionally filters by source cell line(s).
#' 1. Retrieves each source signature and filters genes by direction and magnitude.
#' 1. Queries iLINCS for concordant signatures in the chosen output library.
#' 1. Computes paired or unpaired consensus concordance across up/down regulated sets.
#' 1. Returns an augmented tibble of similarity scores and rich source/target metadata.
#'
#' The paired workflow evaluates concordance separately for up‑ and down‑regulated
#' genes and then combines (via [consensusConcordants()]) the two result sets. When
#' `paired = FALSE` a single aggregate signature (direction = "any") is used.
#'
#' Network access: This function performs remote API calls (unless tests are run
#' under a mocking context such as `httptest2::with_mock_api()`). Examples are
#' wrapped in `\donttest{}` to avoid false negatives on CRAN / Bioconductor builders
#' without network access.
#'
#' @section Thresholds:
#' * `filterThreshold` controls gene selection within each source signature. It is
#'   passed to [filterSignature()] as the absolute (or directional) threshold.
#' * `similarityThreshold` is applied when forming the consensus concordants to
#'   discard low similarity entries.
#'
#' @param target Character scalar. Gene symbol (for KD/OE libraries), or drug /
#'   compound name (for CP library) used to locate source signatures.
#' @param inputLib Character ("OE", "KD", or "CP"). Library from which source
#'   signatures for `target` are drawn.
#' @param outputLib Character ("OE", "KD", or "CP"). Library queried for
#'   concordant signatures.
#' @param filterThreshold Numeric in \(0,1\]. Minimum absolute (or directional)
#'   change used to retain genes in each source signature prior to concordance.
#'   Default `0.85` is conservative; consider lowering (e.g. `0.5`) for broader
#'   coverage.
#' @param similarityThreshold Numeric in \[0,1\]. Minimum similarity score retained
#'   in the final consensus result set. Default `0.321` (~ upper third).
#' @param paired Logical. If `TRUE` (default) computes concordance separately for
#'   up and down regulated gene sets; if `FALSE` uses all selected genes together.
#' @param inputCellLines Optional character vector restricting the search for
#'   source signatures to specified cell line(s). If `NULL` all available are
#'   considered.
#' @param outputCellLines Optional character vector restricting target signatures
#'   (during consensus formation) to specified cell line(s). If `NULL` all are
#'   considered.
#'
#' @return A tibble (data frame) with one row per consensus concordant target
#'   signature. Typical columns include:
#'   * `Source` / `Target` – gene or compound names.
#'   * `Similarity` – numeric concordance score in \[-1,1\].
#'   * `SourceSignature`, `TargetSignature` – iLINCS signature identifiers.
#'   * `SourceCellLine`, `TargetCellLine` – originating cell lines (if applicable).
#'   * `SourceConcentration`, `TargetConcentration` – dosing information for CP.
#'   * `SourceTime`, `TargetTime` – time point metadata.
#'
#' @details Errors are raised if:
#' * No source signatures match `target` in the requested `inputLib` (empty set).
#' * Invalid library codes are supplied.
#'
#' Internally this function orchestrates: [getSignature()], [filterSignature()],
#' [getConcordants()] and [consensusConcordants()]. It returns a vertically
#' concatenated result across all matching source signatures.
#'
#' @seealso [getSignature()], [filterSignature()], [getConcordants()],
#'   [consensusConcordants()], [prepareSignature()] for lower‑level operations.
#'
#' @export
#'
#' @importFrom dplyr filter pull select any_of inner_join mutate
#' @importFrom stringr str_to_lower
#' @importFrom purrr map map2 map_dfr
#' @importFrom rlang .data
#'
#' @examples
#' # Input validation examples (no API calls)
#' # Demonstrate library parameter validation
#' tryCatch(
#'     investigateTarget(target = "TP53", inputLib = "INVALID", outputLib = "CP"),
#'     error = function(e) message("Expected error: invalid inputLib")
#' )
#'
#' tryCatch(
#'     investigateTarget(target = "TP53", inputLib = "KD", outputLib = "INVALID"),
#'     error = function(e) message("Expected error: invalid outputLib")
#' )
#'
#' \donttest{
#' # This function makes multiple API calls to iLINCS and may take several minutes
#' # Basic paired investigation of a knockdown signature against compound library
#' set.seed(1)
#' res <- investigateTarget(
#'     target = "AATK",
#'     inputLib = "KD",
#'     outputLib = "CP",
#'     filterThreshold = 0.5,
#'     similarityThreshold = 0.3,
#'     paired = TRUE
#' )
#' head(res)
#'
#' # Unpaired (aggregate) workflow — often faster, returns a single consensus set
#' res_unpaired <- investigateTarget(
#'     target = "AATK", inputLib = "KD", outputLib = "CP",
#'     filterThreshold = 0.5, similarityThreshold = 0.3, paired = FALSE
#' )
#' head(res_unpaired)
#'
#' # Restrict source signatures to specific cell lines (if available)
#' # and target signatures to a subset of cell lines during consensus
#' res_filtered <- investigateTarget(
#'     target = "AATK", inputLib = "KD", outputLib = "CP",
#'     outputCellLines = c("MCF7"),
#'     filterThreshold = 0.5, similarityThreshold = 0.3
#' )
#' head(res_filtered)
#'
#' # Using httptest2 (if installed) to mock network calls:
#' # httptest2::with_mock_api({
#' #   mock_res <- investigateTarget("AATK", "KD", "CP", filterThreshold = 0.5)
#' #   print(head(mock_res))
#' # })
#' }
investigateTarget <- function(
    target,
    inputLib, outputLib,
    filterThreshold = 0.85,
    similarityThreshold = 0.321,
    paired = TRUE, inputCellLines = NULL,
    outputCellLines = NULL) {
    stopIfInvalidLibraries(c(inputLib, outputLib)) # nolint: object_usage_linter.

    # Load metadata and obtain candidate source signatures for the target
    inputMetadata <- .loadMetadata(inputLib) # nolint: object_usage_linter.
    filteredSignatureIds <- inputMetadata |>
        dplyr::filter(
            stringr::str_to_lower(.data[["Source"]]) == # nolint: object_usage_linter.
                stringr::str_to_lower(target)
        ) |>
        (\(x) if (!is.null(inputCellLines)) dplyr::filter(x, .data[["SourceCellLine"]] %in% inputCellLines) else x)() |>
        dplyr::pull(.data[["SourceSignature"]])

    if (length(filteredSignatureIds) == 0L) {
        stop("No signatures match the given input criteria.", call. = FALSE)
    }

    # Helper to process a single signature ID using existing modular functions
    processOne <- function(sigId) {
        sig <- getSignature(sigId) # nolint: object_usage_linter.
        .computeConsensusFromSignature( # nolint: object_usage_linter.
            sig,
            outputLib = outputLib,
            filterThreshold = filterThreshold,
            similarityThreshold = similarityThreshold,
            paired = paired,
            outputCellLines = outputCellLines
        ) |>
            dplyr::mutate(SourceSignature = sigId)
    }

    # Process all source signatures and augment with metadata
    filteredSignatureIds |>
        purrr::map(processOne) |>
        dplyr::bind_rows() |>
        dplyr::inner_join(inputMetadata, by = "SourceSignature") |>
        dplyr::select(
            dplyr::any_of(c(
                "Source", "Target", "Similarity", "SourceSignature",
                "SourceCellLine", "SourceConcentration", "SourceTime",
                "TargetSignature", "TargetCellLine", "TargetConcentration", "TargetTime",
                "InputSigDirection", "SignatureType", "pValue"
            ))
        )
}
