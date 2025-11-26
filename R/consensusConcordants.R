#' @include utilities.R
NULL

#' Validate consensusConcordants input parameters
#'
#' This internal function validates all input parameters for the consensusConcordants
#' function to ensure they meet the required constraints.
#'
#' @param dots A list of dataframes passed via ... parameter.
#' @param paired Logical indicating whether paired analysis is requested.
#' @param cutoff Numeric similarity cutoff value.
#' @param cellLine Character vector of cell lines, or NULL.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' This function performs the following validations:
#' \enumerate{
#'   \item Ensures paired analysis has exactly two dataframes
#'   \item Ensures unpaired analysis has exactly one dataframe
#'   \item Validates cutoff is numeric and within reasonable range
#'   \item Validates cellLine parameter format
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid calls (no errors)
#' testData <- data.frame(similarity = c(0.5, -0.3), compound = c("A", "B"))
#' .validateConsensusConcordantsInput(list(testData), FALSE, 0.3, NULL)
#' .validateConsensusConcordantsInput(list(testData, testData), TRUE, 0.3, "A375")
#'
#' # Invalid calls (will throw errors)
#' .validateConsensusConcordantsInput(list(), FALSE, 0.3, NULL) # No data
#' .validateConsensusConcordantsInput(list(testData), TRUE, 0.3, NULL) # Paired needs 2 dataframes
#' }
# nolint start: object_length_linter, cyclocomp_linter.
.validateConsensusConcordantsInput <- function(dots, paired, cutoff, cellLine) {
    # nolint end.
    # Validate input count based on analysis type
    if (paired && length(dots) != 2L) {
        stop("Paired analysis requires two data frames", call. = FALSE)
    }

    if (!paired && length(dots) != 1L) {
        stop("Unpaired analysis requires only one dataframe", call. = FALSE)
    }

    # Validate cutoff parameter
    if (!is.numeric(cutoff) || length(cutoff) != 1L) {
        stop("cutoff must be a single numeric value", call. = FALSE)
    }

    if (cutoff < 0L || cutoff > 1L) {
        stop("cutoff must be between 0 and 1", call. = FALSE)
    }

    # Validate cellLine parameter
    if (!is.null(cellLine) && !is.character(cellLine)) {
        stop("cellLine must be a character vector or NULL", call. = FALSE)
    }

    # Validate that dataframes are not empty and have required columns
    for (i in seq_along(dots)) {
        concordantsDf <- dots[[i]]

        if (is.null(concordantsDf) || nrow(concordantsDf) == 0L) {
            stop("All input dataframes must be non-empty", call. = FALSE)
        }

        # Check for essential columns
        requiredCols <- c("treatment", "similarity")
        missingCols <- setdiff(requiredCols, names(concordantsDf))
        if (length(missingCols) > 0L) {
            stop("Missing required columns in dataframe ", i, ": ",
                toString(missingCols),
                call. = FALSE
            )
        }
    }
}

#' Combine concordants dataframes for consensus analysis
#'
#' This internal function combines one or more concordants dataframes into a
#' single dataframe for further processing.
#'
#' @param dots A list of dataframes to combine.
#'
#' @return A combined dataframe with all input data.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Combines multiple dataframes using row binding
#'   \item Preserves all columns from input dataframes
#'   \item Handles cases where dataframes have different column sets
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' df1 <- data.frame(similarity = 0.5, compound = "A")
#' df2 <- data.frame(similarity = -0.3, compound = "B")
#' combined <- .combineConcordantsData(list(df1, df2))
#' }
.combineConcordantsData <- function(dots) {
    dplyr::bind_rows(dots)
}

#' Filter concordants data by cell line
#'
#' This internal function filters the concordants data to include only
#' the specified cell lines.
#'
#' @param concordants A dataframe containing concordants data.
#' @param cellLine A character vector of cell lines to include, or NULL for no filtering.
#'
#' @return A filtered dataframe containing only the specified cell lines.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Filters data based on the cellline column
#'   \item Returns original data if cellLine is NULL
#'   \item Handles cases where no data matches the specified cell lines
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' testData <- data.frame(
#'     similarity = c(0.5, -0.3, 0.7),
#'     cellline = c("A375", "PC3", "MCF7")
#' )
#' filtered <- .filterByCellLine(testData, c("A375", "PC3"))
#' }
.filterByCellLine <- function(concordants, cellLine) {
    if (is.null(cellLine)) {
        return(concordants)
    }

    dplyr::filter(concordants, .data[["cellline"]] %in% cellLine) # nolint: object_usage_linter.
}

#' Apply similarity cutoff filter to concordants data
#'
#' This internal function filters concordants data based on absolute similarity values
#' meeting or exceeding the specified cutoff threshold.
#'
#' @param concordants A dataframe containing concordants data.
#' @param cutoff Numeric similarity cutoff value.
#'
#' @return A filtered dataframe containing only entries meeting the similarity cutoff.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Filters based on absolute similarity values
#'   \item Retains both positive and negative similarities above threshold
#'   \item Removes entries below the cutoff threshold
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' testData <- data.frame(similarity = c(0.5, -0.8, 0.2, -0.1))
#' filtered <- .applySimilarityCutoff(testData, 0.3)
#' # Returns entries with |similarity| >= 0.3
#' }
.applySimilarityCutoff <- function(concordants, cutoff) {
    dplyr::filter(concordants, abs(.data[["similarity"]]) >= cutoff) # nolint: object_usage_linter.
}

#' Group concordants by target and select maximum similarity entries
#'
#' This internal function groups concordants data by target (compound or treatment)
#' and retains only the entries with maximum absolute similarity for each target.
#'
#' @param concordants A dataframe containing filtered concordants data.
#'
#' @return A dataframe with deduplicated targets, keeping maximum similarity entries.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Groups by treatment or compound columns (whichever is available)
#'   \item For each group, retains only entries with maximum absolute similarity
#'   \item Handles ties by keeping all tied entries
#'   \item Preserves the structure for downstream processing
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr group_by across any_of filter ungroup
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' testData <- data.frame(
#'     compound = c("A", "A", "B", "B"),
#'     similarity = c(0.5, 0.8, -0.3, -0.7),
#'     cellline = c("A375", "PC3", "A375", "PC3")
#' )
#' grouped <- .groupByTargetAndSelectMax(testData)
#' # Returns entries with max |similarity| for each compound
#' }
.groupByTargetAndSelectMax <- function(concordants) {
    if (nrow(concordants) == 0L) {
        return(concordants)
    }

    concordants |>
        dplyr::group_by(
            dplyr::across(!!"treatment")
        ) |>
        dplyr::filter(
            abs(.data[["similarity"]]) == max(abs(.data[["similarity"]])) # nolint: object_usage_linter.
        ) |>
        dplyr::ungroup()
}

#' Select and order consensus results columns
#'
#' This internal function selects the relevant columns for consensus results
#' and orders them appropriately for output.
#'
#' @param concordants A dataframe containing processed concordants data.
#'
#' @return A dataframe with selected and ordered columns for consensus output.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Selects standard consensus output columns
#'   \item Orders results by descending absolute similarity
#'   \item Handles both CP/KD libraries (with concentration) and OE libraries (without)
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr select any_of arrange
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' testData <- data.frame(
#'     signatureid = "SIG1",
#'     compound = "A",
#'     cellline = "A375",
#'     similarity = 0.8,
#'     sig_direction = "Up",
#'     pValue = 0.01
#' )
#' selected <- .selectAndOrderResults(testData)
#' }
.selectAndOrderResults <- function(concordants) {
    concordants |>
        dplyr::select(
            dplyr::any_of(c(
                "signatureid", "treatment", "cellline", "time",
                "concentration", "sig_direction", "sig_type",
                "similarity", "pValue"
            ))
        ) |>
        dplyr::arrange(dplyr::desc(abs(.data[["similarity"]]))) # nolint: object_usage_linter.
}

#' Apply target column renaming to consensus results
#'
#' This internal function applies the standard target column renaming
#' to produce the final consensus concordants output format.
#'
#' @param concordants A dataframe containing selected consensus results.
#'
#' @return A dataframe with renamed columns following consensus output standards.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Applies targetRename function to standardize column names
#'   \item Converts internal column names to user-facing consensus format
#'   \item Handles different library types appropriately
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr rename_with
#'
#' @examples
#' \dontrun{
#' testData <- data.frame(
#'     signatureid = "SIG1",
#'     compound = "A",
#'     cellline = "A375",
#'     similarity = 0.8
#' )
#' renamed <- .applyTargetRenaming(testData)
#' }
.applyTargetRenaming <- function(concordants) {
    dplyr::rename_with(concordants, targetRename) # nolint: object_usage_linter.
}

#' Process concordants data through the complete consensus pipeline
#'
#' This internal function orchestrates the complete processing pipeline for
#' consensus concordants analysis.
#'
#' @param concordants A combined dataframe containing all concordants data.
#' @param cutoff Numeric similarity cutoff value.
#' @param cellLine Character vector of cell lines to include, or NULL.
#'
#' @return A processed dataframe with consensus concordants results.
#'
#' @details
#' This function coordinates the following processing steps:
#' \enumerate{
#'   \item Cell line filtering (if specified)
#'   \item Similarity cutoff application
#'   \item Target grouping and maximum similarity selection
#'   \item Column selection and ordering
#'   \item Target column renaming
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' testData <- data.frame(
#'     similarity = c(0.5, -0.8, 0.2),
#'     compound = c("A", "B", "C"),
#'     cellline = c("A375", "PC3", "A375")
#' )
#' processed <- .processConsensusPipeline(testData, 0.3, "A375")
#' }
.processConsensusPipeline <- function(concordants, cutoff, cellLine) {
    concordants |>
        .filterByCellLine(cellLine) |>
        .applySimilarityCutoff(cutoff) |>
        .groupByTargetAndSelectMax() |>
        .selectAndOrderResults() |>
        .applyTargetRenaming()
}

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
#' @importFrom dplyr filter arrange any_of group_by across
#' @importFrom dplyr select bind_rows rename_with ungroup
#' @importFrom rlang .data
#'
#' @examples
#' # Create mock concordants data for demonstration
#' mockConcordants <- data.frame(
#'     signatureid = paste0("SIG", 1:10),
#'     treatment = c(
#'         "TP53", "TP53", "MYC", "MYC", "EGFR",
#'         "EGFR", "KRAS", "BRCA1", "BRCA1", "PIK3CA"
#'     ),
#'     cellline = c(
#'         "A375", "PC3", "A375", "MCF7", "A375",
#'         "PC3", "A375", "A375", "MCF7", "A375"
#'     ),
#'     time = rep("24H", 10),
#'     concentration = rep(NA, 10),
#'     sig_direction = rep("DOWN", 10),
#'     sig_type = rep("single", 10),
#'     similarity = c(
#'         0.85, 0.72, -0.68, -0.45, 0.55,
#'         0.38, 0.42, 0.51, 0.33, 0.29
#'     ),
#'     pValue = rep(0.001, 10)
#' )
#'
#' # Example 1: Basic consensus with default cutoff
#' consensus <- consensusConcordants(mockConcordants)
#' nrow(consensus) # Targets with |similarity| >= 0.321
#'
#' # Example 2: Consensus with higher cutoff
#' consensus_strict <- consensusConcordants(mockConcordants, cutoff = 0.5)
#' nrow(consensus_strict) # Fewer targets with higher threshold
#'
#' # Example 3: Filter by cell line
#' consensus_A375 <- consensusConcordants(mockConcordants, cellLine = "A375")
#' unique(consensus_A375$CellLine) # Only A375
#'
#' \donttest{
#' # Network-dependent examples using real iLINCS data
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
#'
#' # Get concordant gene knockdown signatures
#' concordantSignatures <- getConcordants(kdSignature, ilincsLibrary = "KD")
#'
#' # Get the consensus list with different parameters
#' consensus <- consensusConcordants(concordantSignatures, cutoff = 0.5)
#'
#' # Paired analysis example
#' filteredUp <- filterSignature(kdSignature, direction = "up", threshold = 0.5)
#' filteredDown <- filterSignature(kdSignature, direction = "down", threshold = -0.5)
#' concordants_up <- getConcordants(filteredUp, ilincsLibrary = "KD")
#' concordants_down <- getConcordants(filteredDown, ilincsLibrary = "KD")
#' consensus <- consensusConcordants(concordants_up, concordants_down, paired = TRUE)
#' }
#'
consensusConcordants <- function(
  ...,
  paired = FALSE,
  cutoff = 0.321,
  cellLine = NULL
) {
    # Capture input dataframes
    dots <- list(...)

    # Validate all input parameters
    .validateConsensusConcordantsInput(dots, paired, cutoff, cellLine)

    # Combine input dataframes and process through consensus pipeline
    .combineConcordantsData(dots) |>
        .processConsensusPipeline(cutoff, cellLine)
}
