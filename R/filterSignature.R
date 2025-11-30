#' Validate filterSignature input parameters
#'
#' This internal function validates all input parameters for the filterSignature
#' function to ensure they meet the required constraints and are mutually compatible.
#'
#' @param signature A data.frame-like object (data.frame, tibble, or DataFrame)
#'   containing the L1000 signature data.
#' @param direction Character string specifying the filtering direction.
#'   Must be one of "up", "down", or "any".
#' @param threshold Numeric value or vector specifying absolute threshold(s).
#'   Can be NULL, a single value, or a vector of two values. Cannot be specified
#'   together with \code{prop}.
#' @param prop Numeric value specifying the proportion for quantile-based filtering.
#'   Must be between 0 and 1. Cannot be specified together with \code{threshold}.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' This function performs the following validations in order:
#' \enumerate{
#'   \item Ensures \code{signature} is a data.frame-like object
#'   \item Validates \code{direction} is one of the allowed values
#'   \item Verifies that only one of \code{threshold} or \code{prop} is specified
#'   \item For \code{threshold}: checks length (1-2 values) and order (lower, higher)
#'   \item For \code{prop}: checks it's a single value between 0 and 1
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid calls (no errors)
#' sig <- data.frame(Value_LogDiffExp = c(-2, -1, 0, 1, 2))
#' .validateFilterSignatureInput(sig, "any", 1.0, NULL)
#' .validateFilterSignatureInput(sig, "up", NULL, 0.1)
#' .validateFilterSignatureInput(sig, "down", c(-1.5, 1.0), NULL)
#'
#' # Invalid calls (will throw errors)
#' .validateFilterSignatureInput(sig, "invalid", 1.0, NULL) # Invalid direction
#' .validateFilterSignatureInput(sig, "any", 1.0, 0.1) # Both threshold and prop
#' .validateFilterSignatureInput(sig, "any", NULL, NULL) # Neither threshold nor prop
#' .validateFilterSignatureInput(sig, "any", c(1, 2, 3), NULL) # Too many thresholds
#' .validateFilterSignatureInput(sig, "any", c(2, 1), NULL) # Wrong threshold order
#' .validateFilterSignatureInput(sig, "any", NULL, 1.5) # Proportion > 1
#' .validateFilterSignatureInput(sig, "any", NULL, -0.1) # Proportion < 0
#' }
.validateFilterSignatureInput <- function(signature, direction, threshold, prop) { # nolint: cyclocomp_linter.
    # 1. Validate signature data structure
    if (!inherits(signature, c("data.frame", "DFrame"))) {
        stop("Signature must be a data.frame, tibble, or DataFrame", call. = FALSE)
    }

    # 2. Validate direction parameter
    if (!direction %in% c("up", "down", "any")) {
        stop("Direction must be one of 'up', 'down' or 'any'", call. = FALSE)
    }

    # 3. Validate threshold/prop mutual exclusivity
    if (!is.null(threshold) && !is.null(prop)) {
        stop("Only one of prop or threshold can be specified", call. = FALSE)
    }
    if (is.null(threshold) && is.null(prop)) {
        stop("One of prop or threshold must be specified", call. = FALSE)
    }

    # 4. Validate threshold parameter
    if (!is.null(threshold)) {
        if (length(threshold) < 1L || length(threshold) > 2L) {
            stop("Threshold must be specified as one or two values", call. = FALSE)
        }
        if (length(threshold) == 2L && threshold[1L] > threshold[2L]) {
            stop(
                "When two thresholds are specified, they must be in order (lower, higher)",
                call. = FALSE
            )
        }
    }

    # 5. Validate proportion parameter
    if (!is.null(prop)) {
        if (length(prop) != 1L) {
            stop("Proportion must be specified as a single value", call. = FALSE)
        }
        if (prop <= 0L) {
            stop("Proportion must be between greater than 0", call. = FALSE)
        }
        if (prop > 0.5) {
            stop("Proportion must be between less than 0.5", call. = FALSE)
        }
    }
}

#' Calculate thresholds from single threshold value
#'
#' This internal function creates symmetric filtering thresholds from a single
#' threshold value. The input value is used as the positive threshold, and its
#' negative is used as the negative threshold.
#'
#' @param threshold A single positive numeric value representing the absolute
#'   threshold for filtering.
#'
#' @return A named list with two elements:
#'  \itemize{
#'   \item `downThreshold`: The negative threshold (-threshold)
#'   \item `upThreshold`: The positive threshold (threshold)
#' }
#'
#' @details
#' This function is used when a single threshold value is provided to
#' \code{filterSignature}. It creates symmetric thresholds where genes with
#' log fold-change values greater than or equal to the positive threshold
#' (up-regulated) or less than or equal to the negative threshold
#' (down-regulated) are retained.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Create symmetric thresholds from threshold = 1.5
#' thresholds <- .calculateSingleThreshold(1.5)
#' # Returns: list(downThreshold = -1.5, upThreshold = 1.5)
#'
#' # Create symmetric thresholds from threshold = 0.8
#' thresholds <- .calculateSingleThreshold(0.8)
#' # Returns: list(downThreshold = -0.8, upThreshold = 0.8)
#' }
.calculateSingleThreshold <- function(threshold) {
    list(
        downThreshold = -threshold,
        upThreshold = threshold
    )
}

#' Calculate thresholds from two threshold values
#'
#' This internal function handles asymmetric filtering thresholds when two
#' threshold values are provided. The first value is used as the down-regulated
#' threshold and the second value is used as the up-regulated threshold.
#'
#' @param threshold A numeric vector of length 2 containing the threshold values.
#'   The first element is the down-regulated threshold (typically negative),
#'   and the second element is the up-regulated threshold (typically positive).
#'
#' @return A named list with two elements:
#'  \itemize{
#'   \item `downThreshold`: The down-regulated threshold (`threshold[1]`)
#'   \item `upThreshold`: The up-regulated threshold (`threshold[2]`)
#' }
#'
#' @details
#' This function enables asymmetric filtering where different absolute thresholds
#' can be applied to up-regulated and down-regulated genes. This is useful when
#' you want to apply stricter criteria to one direction of regulation than the other.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Create asymmetric thresholds
#' thresholds <- .calculateDoubleThreshold(c(-2.0, 1.5))
#' # Returns: list(downThreshold = -2.0, upThreshold = 1.5)
#'
#' # Stricter threshold for down-regulation
#' thresholds <- .calculateDoubleThreshold(c(-1.0, 0.5))
#' # Returns: list(downThreshold = -1.0, upThreshold = 0.5)
#'
#' # Equal but explicit thresholds
#' thresholds <- .calculateDoubleThreshold(c(-1.5, 1.5))
#' # Returns: list(downThreshold = -1.5, upThreshold = 1.5)
#' }
.calculateDoubleThreshold <- function(threshold) {
    if (threshold[[1L]] > threshold[[2L]]) {
        stop(
            "When two thresholds are specified, they must be in order (lower, higher)",
            call. = FALSE
        )
    }
    list(
        downThreshold = threshold[[1L]],
        upThreshold = threshold[[2L]]
    )
}

#' Calculate thresholds using absolute threshold values
#'
#' This internal function coordinates the calculation of filtering thresholds
#' when absolute threshold values are provided. It dispatches to the appropriate
#' calculation function based on the number of threshold values provided.
#'
#' @param threshold A numeric value or vector specifying the absolute threshold(s).
#'
#'   Can be:
#'
#'     - A single value: Dispatched to [.calculateSingleThreshold()]
#'     - A vector of two values: Dispatched to [.calculateDoubleThreshold()]
#'
#' @return A named list with two elements:
#'   - `downThreshold`: The threshold for down-regulated genes
#'   - `upThreshold`: The threshold for up-regulated genes
#'
#' @details
#' This function serves as a dispatcher that:
#'   - Checks the length of the threshold parameter
#'   - Calls the appropriate threshold calculation function
#'   - Throws an error if an invalid number of thresholds is provided
#'
#' The function ensures that only single values or pairs of values are accepted,
#' maintaining the integrity of the filtering logic.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Single threshold - creates symmetric thresholds
#' thresholds <- .calculateAbsoluteThresholds(1.0)
#' # Returns: list(downThreshold = -1.0, upThreshold = 1.0)
#'
#' # Double threshold - uses provided values
#' thresholds <- .calculateAbsoluteThresholds(c(-1.5, 2.0))
#' # Returns: list(downThreshold = -1.5, upThreshold = 2.0)
#'
#' # Invalid - too many values (will throw error)
#' # thresholds <- .calculateAbsoluteThresholds(c(1.0, 2.0, 3.0))
#' }
.calculateAbsoluteThresholds <- function(threshold) {
    if (length(threshold) == 1L) {
        .calculateSingleThreshold(threshold)
    } else if (length(threshold) == 2L) {
        .calculateDoubleThreshold(threshold)
    } else {
        stop("Threshold must be specified as one or two values", call. = FALSE)
    }
}

#' Calculate thresholds using proportional values
#'
#' This internal function calculates filtering thresholds based on quantiles
#' of the log fold-change distribution in the signature data. This enables
#' proportion-based filtering that adapts to the data distribution.
#'
#' @param signature A data.frame-like object containing the signature data.
#'   Must have a column named "Value_LogDiffExp" containing log fold-change values.
#' @param prop A numeric value between 0 and 1 specifying the proportion of
#'   genes to select from each tail of the distribution.
#'
#' @return A named list with two elements:
#'   \item{downThreshold}{The quantile threshold for down-regulated genes (quantile at prop)}
#'   \item{upThreshold}{The quantile threshold for up-regulated genes (quantile at 1-prop)}
#'
#' @details
#' This function calculates thresholds using the \code{quantile} function:
#' \itemize{
#'   \item \code{downThreshold}: The \code{prop} quantile of the expression values
#'   \item \code{upThreshold}: The \code{1-prop} quantile of the expression values
#' }
#'
#' For example, with \code{prop = 0.1}:
#' \itemize{
#'   \item \code{downThreshold}: 10th percentile (bottom 10% of values)
#'   \item \code{upThreshold}: 90th percentile (top 10% of values)
#' }
#'
#' This approach is particularly useful when you want to select a fixed proportion
#' of the most differentially expressed genes regardless of their absolute
#' expression values.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Create sample signature data
#' signature <- data.frame(
#'     Value_LogDiffExp = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6)
#' )
#'
#' # Calculate thresholds for top/bottom 20%
#' thresholds <- .calculateProportionalThreshold(signature, 0.2)
#' # Returns thresholds based on 20th and 80th percentiles
#'
#' # Calculate thresholds for top/bottom 10%
#' thresholds <- .calculateProportionalThreshold(signature, 0.1)
#' # Returns thresholds based on 10th and 90th percentiles
#'
#' # Calculate thresholds for top/bottom 5% (most extreme)
#' thresholds <- .calculateProportionalThreshold(signature, 0.05)
#' # Returns thresholds based on 5th and 95th percentiles
#' }
.calculateProportionalThreshold <- function(signature, prop) { # nolint: object_length_linter.
    limits <- round(
        quantile(
            signature[["Value_LogDiffExp"]], c(prop, 1L - prop)
        ), 2L
    )

    list(
        downThreshold = limits[1L],
        upThreshold = limits[2L]
    )
}

#' Apply filtering based on direction and thresholds
#'
#' This internal function performs the actual filtering of the signature data
#' based on the specified direction and calculated thresholds. It implements
#' the core filtering logic using dplyr operations.
#'
#' @param signature A data.frame-like object containing the signature data.
#'   Must have a column named "Value_LogDiffExp" containing log fold-change values.
#' @param direction Character string specifying the filtering direction.
#'   Must be one of:
#'   \itemize{
#'     \item "up": Keep only up-regulated genes (logFC >= upThreshold)
#'     \item "down": Keep only down-regulated genes (logFC <= downThreshold)
#'     \item "any": Keep both up- and down-regulated genes (logFC >= upThreshold OR logFC <= downThreshold)
#'   }
#' @param thresholds A named list containing:
#'   \itemize{
#'   \item `downThreshold`: Threshold for down-regulated genes
#'   \item `upThreshold`: Threshold for up-regulated genes
#' }
#' @return A tibble containing the filtered signature data with the same structure
#'   as the input but including only rows that meet the filtering criteria.
#'
#' @details
#' The filtering logic depends on the direction parameter:
#' \itemize{
#'   \item \strong{"up"}: Retains genes where \code{Value_LogDiffExp >= upThreshold}
#'   \item \strong{"down"}: Retains genes where \code{Value_LogDiffExp <= downThreshold}
#'   \item \strong{"any"}: Retains genes where
#'     \code{Value_LogDiffExp >= upThreshold OR Value_LogDiffExp <= downThreshold}
#' }
#'
#' The function uses \code{dplyr::filter} with \code{rlang::.data} for
#' non-standard evaluation, ensuring compatibility with different data frame types
#' and avoiding issues with variable scoping.
#'
#' @keywords internal
#'
#' @import DFplyr
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom stats quantile
#'
#'
#' @examples
#' \dontrun{
#' # Create sample signature data
#' signature <- data.frame(
#'     signatureID = rep("TEST", 10),
#'     Name_GeneSymbol = paste0("GENE", 1:10),
#'     Value_LogDiffExp = c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4)
#' )
#'
#' # Define thresholds
#' thresholds <- list(downThreshold = -1.5, upThreshold = 1.5)
#'
#' # Filter for up-regulated genes only
#' up_filtered <- .applyDirectionFilter(signature, "up", thresholds)
#' # Returns genes with logFC >= 1.5 (GENE8, GENE9, GENE10)
#'
#' # Filter for down-regulated genes only
#' down_filtered <- .applyDirectionFilter(signature, "down", thresholds)
#' # Returns genes with logFC <= -1.5 (GENE1, GENE2)
#'
#' # Filter for both up- and down-regulated genes
#' both_filtered <- .applyDirectionFilter(signature, "any", thresholds)
#' # Returns genes with |logFC| >= 1.5 (GENE1, GENE2, GENE8, GENE9, GENE10)
#' }
.applyDirectionFilter <- function(signature, thresholds, direction = "any") {
    downThreshold <- thresholds[["downThreshold"]]
    upThreshold <- thresholds[["upThreshold"]]

    if (direction == "up") {
        dplyr::filter(signature, .data[["Value_LogDiffExp"]] >= upThreshold) # nolint: object_usage_linter.
    } else if (direction == "down") {
        dplyr::filter(signature, .data[["Value_LogDiffExp"]] <= downThreshold) # nolint: object_usage_linter.
    } else {
        dplyr::filter(
            signature,
            .data[["Value_LogDiffExp"]] >= !!upThreshold | # nolint: object_usage_linter.
                .data[["Value_LogDiffExp"]] <= !!downThreshold
        )
    }
}

#' Filter the L1000 Signature
#' `r lifecycle::badge("stable")`
#'
#' This function filters the L1000 signature to a given threshold, identifying
#' up-regulated, down-regulated, or both up- and down-regulated genes. The
#' function supports both absolute threshold filtering and proportional filtering
#' based on quantiles of the expression data.
#'
#' @param signature A data.frame, tibble, or DataFrame containing the L1000 signature.
#'   Must contain a column named "Value_LogDiffExp" with log fold-change values.
#' @param direction Character string specifying the direction to filter.
#'   Must be one of "up" (up-regulated genes only), "down" (down-regulated genes only),
#'   or "any" (both up- and down-regulated genes). Defaults to "any".
#' @param threshold Numeric value or vector specifying the log fold-change threshold(s).
#'   Can be:
#'   \itemize{
#'     \item A single positive value: Creates symmetric thresholds (\eqn{\pm threshold})
#'     \item A vector of two values: First value is the down-regulated threshold,
#'           second value is the up-regulated threshold
#'   }
#'   Cannot be specified together with \code{prop}. One of \code{threshold} or
#'   \code{prop} must be provided.
#' @param prop Numeric value between 0 and 1 specifying the proportion of genes
#'   to select from the top and bottom of the expression distribution. For example,
#'   \code{prop = 0.1} selects the top 10% most up-regulated and bottom 10%
#'   most down-regulated genes. Cannot be specified together with \code{threshold}.
#'
#' @return A tibble containing the filtered L1000 signature with the same structure
#'   as the input but containing only genes that meet the filtering criteria.
#'
#' @details
#' The filtering process follows these steps:
#' \enumerate{
#'   \item Input validation: Checks data frame structure and parameter consistency
#'   \item Threshold calculation: Computes filtering thresholds based on either
#'         absolute values (\code{threshold}) or quantiles (\code{prop})
#'   \item Direction-based filtering: Applies the computed thresholds according
#'         to the specified direction
#' }
#'
#' When using \code{threshold}:
#' \itemize{
#'   \item Single value: Genes with |logFC| >= threshold are retained
#'   \item Two values: Genes with logFC <= `threshold[1]` OR logFC >= `threshold[2]`
#' }
#'
#' When using \code{prop}:
#' \itemize{
#'   \item Thresholds are calculated as quantiles of the expression distribution
#'   \item Down threshold = quantile(logFC, prop)
#'   \item Up threshold = quantile(logFC, 1 - prop)
#' }
#'
#' @seealso
#' \code{\link{getSignature}} for retrieving L1000 signatures from iLINCS,
#' \code{\link{prepareSignature}} for preparing custom signatures,
#' \code{\link{getConcordants}} for finding concordant signatures
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @import DFplyr
#'
#' @examples
#' # Create a mock signature for demonstration
#' mockSignature <- data.frame(
#'     signatureID = rep("MOCK001", 20),
#'     Name_GeneSymbol = paste0("GENE", 1:20),
#'     ID_geneid = 1:20,
#'     Value_LogDiffExp = c(
#'         -3.5, -2.8, -2.1, -1.5, -1.2, -0.8, -0.5, -0.3,
#'         -0.1, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6, 2.0, 2.4, 2.9, 3.3, 3.8
#'     )
#' )
#'
#' # Example 1: Filter by symmetric absolute threshold
#' # Keeps genes with |logFC| >= 1.5
#' filteredSymmetric <- filterSignature(mockSignature, threshold = 1.5)
#' nrow(filteredSymmetric) # Should return 8 genes
#'
#' # Example 2: Filter by asymmetric absolute thresholds
#' # Keeps genes with logFC <= -2.0 OR logFC >= 2.5
#' filteredAsymmetric <- filterSignature(mockSignature, threshold = c(-2.0, 2.5))
#' nrow(filteredAsymmetric) # Should return 5 genes
#'
#' # Example 3: Filter by proportion (top and bottom 20%)
#' filteredProportion <- filterSignature(mockSignature, prop = 0.2)
#' nrow(filteredProportion) # Should return 8 genes (4 up + 4 down)
#'
#' # Example 4: Filter only up-regulated genes by threshold
#' upRegulated <- filterSignature(mockSignature, direction = "up", threshold = 1.0)
#' all(upRegulated$Value_LogDiffExp >= 1.0) # Should be TRUE
#'
#' # Example 5: Filter only down-regulated genes by threshold
#' downRegulated <- filterSignature(mockSignature, direction = "down", threshold = 1.0)
#' all(downRegulated$Value_LogDiffExp <= -1.0) # Should be TRUE
#'
#' \donttest{
#' # Network-dependent examples using real iLINCS data
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
#'
#' # Filter for top 5% most extreme genes
#' topExtreme <- filterSignature(kdSignature, prop = 0.05)
#'
#' # Get top 20% most up-regulated genes
#' topUpregulated <- filterSignature(kdSignature, direction = "up", prop = 0.2)
#' }
filterSignature <- function(
  signature, direction = "any",
  threshold = NULL, prop = NULL
) {
    # Validate input parameters
    .validateFilterSignatureInput(signature, direction, threshold, prop)

    # Calculate thresholds based on input type
    if (!is.null(threshold)) {
        thresholds <- .calculateAbsoluteThresholds(threshold)
    } else {
        thresholds <- .calculateProportionalThreshold(signature, prop)
    }

    # Apply filtering based on direction
    .applyDirectionFilter(signature, thresholds, direction)
}
