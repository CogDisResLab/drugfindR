#' Validate getSignature input parameters
#'
#' This internal function validates all input parameters for the getSignature
#' function to ensure they meet the required constraints.
#'
#' @param sigId A character string containing the iLINCS signature ID to retrieve.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' This function performs the following validations:
#' \enumerate{
#'   \item Ensures \code{sigId} is a character vector of length 1
#'   \item Ensures \code{sigId} is not empty or whitespace-only
#'   \item Validates that the signature exists in the metadata tables
#' }
#'
#' @keywords internal
#'
#' @examples NULL
.validateGetSignatureInput <- function(sigId) {
    # 1. Validate sigId parameter
    if (!is.character(sigId)) {
        stop("sigId must be a character string", call. = FALSE)
    }

    if (length(sigId) != 1L) {
        stop("sigId must be a single character string", call. = FALSE)
    }

    if (!nzchar(trimws(sigId))) {
        stop("sigId cannot be empty or consist only of whitespace", call. = FALSE)
    }

    # 2. Validate signature exists in metadata
    if (!.isValidSignatureId(sigId)) {
        stop(
            "Signature ID '", sigId, "' not found in metadata tables.\n",
            "Please verify the signature ID exists in the iLINCS database.",
            call. = FALSE
        )
    }
}

#' Check if a signature ID exists in the metadata tables
#'
#' This internal function validates whether a signature ID exists in any of
#' the metadata tables (CP, KD, or OE).
#'
#' @param sigId A character string or vector containing the signature ID(s) to validate.
#'
#' @return A logical value or vector: TRUE if the signature exists, FALSE otherwise.
#'   Returns a vector of the same length as the input when given a vector.
#'
#' @details
#' This function searches all three metadata tables:
#' \itemize{
#'   \item Chemical Perturbagen (CP) metadata
#'   \item Knockdown (KD) metadata
#'   \item Overexpression (OE) metadata
#' }
#'
#' The function checks the "SourceSignature" column in each metadata table
#' for the provided signature ID(s).
#'
#' @keywords internal
#'
#' @examples NULL
.isValidSignatureId <- function(sigId) {
    # Check all metadata tables for the signature ID
    cpExists <- sigId %in% cpMetadata[["SourceSignature"]] # nolint: object_usage_linter.
    kdExists <- sigId %in% kdMetadata[["SourceSignature"]] # nolint: object_usage_linter.
    oeExists <- sigId %in% oeMetadata[["SourceSignature"]] # nolint: object_usage_linter.

    # Return TRUE if found in any metadata table
    cpExists | kdExists | oeExists
}

#' Create HTTP request for iLINCS signature retrieval
#'
#' This internal function constructs and configures the HTTP request object
#' for retrieving signature data from the iLINCS API.
#'
#' @param sigId A character string containing the iLINCS signature ID to retrieve.
#'
#' @return An httr2 request object configured for the iLINCS downloadSignature endpoint.
#'
#' @details
#' This function builds a complete HTTP request by:
#' \enumerate{
#'   \item Setting the base URL using \code{.ilincsBaseUrl()}
#'   \item Appending the API path: "ilincsR/downloadSignature"
#'   \item Adding query parameters: sigID and noOfTopGenes (set to Inf for all genes)
#'   \item Setting the HTTP method to POST
#'   \item Adding a user agent string using \code{.returnUserAgent()}
#' }
#'
#' The request is configured but not executed - it must be performed using
#' the request execution function.
#'
#' @keywords internal
#'
#' @importFrom httr2 request req_url_path_append req_url_query req_method req_user_agent
#'
#' @examples NULL
.createSignatureRequest <- function(sigId) {
    httr2::request(.ilincsBaseUrl()) |> # nolint: object_usage_linter.
        httr2::req_url_path_append("ilincsR") |>
        httr2::req_url_path_append("downloadSignature") |>
        httr2::req_url_query(sigID = sigId, noOfTopGenes = Inf) |>
        httr2::req_method("GET") |>
        httr2::req_user_agent(.returnUserAgent()) # nolint: object_usage_linter.
}

#' Execute iLINCS API request with error handling
#'
#' This internal function safely executes an httr2 request and captures errors
#' instead of raising them, allowing downstream functions to handle them appropriately.
#'
#' @param request An httr2_request object to be executed.
#' @param verbose Logical indicating whether to display request details.
#'   Default is FALSE.
#'
#' @return An httr2_response object, including error responses that would
#'   normally cause httr2 to raise an error.
#'
#' @details
#' This function configures the request to not raise errors automatically
#' on HTTP error status codes (4xx, 5xx) by using \code{httr2::req_error()}.
#' Instead, error responses are returned as response objects that can be
#' processed by \code{.processSignatureResponse()} to generate appropriate
#' error messages with context.
#'
#' The function handles:
#' \itemize{
#'   \item Network connection errors
#'   \item HTTP error status codes (400, 401, 403, 404, 500, etc.)
#'   \item Timeout errors
#'   \item Other httr2 request failures
#' }
#'
#' @importFrom httr2 req_error req_perform req_options
#'
#' @keywords internal
#'
#' @examples NULL
.executeSignatureRequest <- function(request, verbose = FALSE) {
    # Configure request to not raise errors on HTTP error status codes
    request <- httr2::req_error(request, is_error = function(resp) FALSE)

    # Optionally add verbosity for debugging
    if (verbose) {
        request <- httr2::req_options(request, verbose = TRUE)
    }

    # Execute the request - errors are captured in response object
    httr2::req_perform(request)
}

#' Process successful API response into signature data frame
#'
#' This internal function processes a successful HTTP response from the iLINCS
#' API and converts it into a standardized signature data frame.
#'
#' @param response An httr2 response object from a successful iLINCS API call.
#'
#' @return A tibble containing the signature data with standardized columns:
#'   \itemize{
#'     \item signatureID: The signature identifier
#'     \item ID_geneid: Character gene identifiers
#'     \item Name_GeneSymbol: Gene symbols
#'     \item Value_LogDiffExp: Log fold-change values (rounded to 12 decimal places)
#'     \item Significance_pvalue: P-values (rounded to 12 decimal places)
#'   }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Extracts JSON data from the response body
#'   \item Maps the "signature" elements from the response
#'   \item Flattens the nested structure into a data frame
#'   \item Removes the "PROBE" column (not needed for analysis)
#'   \item Converts gene IDs to character format
#'   \item Rounds numeric values to 12 decimal places for consistency
#'   \item Adds signature metadata including L1000 status
#' }
#'
#' The rounding ensures consistent precision across different platforms and
#' prevents floating-point precision issues in downstream analyses.
#'
#' @keywords internal
#'
#' @importFrom httr2 resp_body_json
#' @importFrom purrr map flatten_dfr
#' @importFrom dplyr select mutate
#' @importFrom rlang .data
#'
#' @examples NULL
.processSuccessfulResponse <- function(response) {
    httr2::resp_body_json(response) |>
        purrr::map("signature") |>
        purrr::flatten_dfr() |>
        dplyr::select(-"PROBE") |>
        dplyr::mutate(
            ID_geneid = as.character(.data[["ID_geneid"]]), # nolint: object_usage_linter.
            Value_LogDiffExp = round(.data[["Value_LogDiffExp"]], 12L),
            Significance_pvalue = round(.data[["Significance_pvalue"]], 12L)
        ) |>
        dplyr::select(
            "signatureID", "ID_geneid", "Name_GeneSymbol",
            "Value_LogDiffExp", "Significance_pvalue"
        )
}

#' Handle API error responses for signature retrieval
#'
#' This internal function processes error responses from the iLINCS API
#' and generates appropriate error messages for signature retrieval failures.
#'
#' @param response An httr2 response object from the iLINCS API.
#'
#' @return This function always stops execution with an error message.
#'
#' @keywords internal
.processSignatureResponseError <- function(response) {
    status <- httr2::resp_status(response)

    # Extract error details from response body
    errorBody <- tryCatch(
        httr2::resp_body_string(response),
        error = function(e) "Unable to retrieve error details"
    )

    # Create error message based on status code
    if (status == 404L) {
        errorMsg <- paste0(
            "iLINCS API request failed (Status 404): Signature Not Found\n",
            "The specified signature ID was not found in the iLINCS database.\n",
            "Please verify the signature ID is correct and try again.\n",
            "Response: ", errorBody
        )
    } else if (status == 400L) {
        errorMsg <- paste0(
            "iLINCS API request failed (Status 400): Bad Request\n",
            "This typically indicates invalid parameters or signature ID format.\n",
            "Please check the signature ID format and try again.\n",
            "Response: ", errorBody
        )
    } else if (status == 500L) {
        errorMsg <- paste0(
            "iLINCS API request failed (Status 500): Internal Server Error\n",
            "This may indicate server issues or data retrieval problems.\n",
            "Please try again later.\n",
            "Response: ", errorBody
        )
    } else {
        errorMsg <- paste0(
            "iLINCS API request failed with status ", status, ".\n",
            "Response: ", errorBody
        )
    }

    stop(errorMsg, call. = FALSE)
}

#' Process iLINCS API response for signature retrieval
#'
#' This internal function dispatches to appropriate handlers based on the
#' response status from the iLINCS API.
#'
#' @param response An httr2 response object from the iLINCS API.
#'
#' @return A tibble containing signature data with standardized columns.
#'
#' @details
#' The function dispatches to specialized handlers:
#' \enumerate{
#'   \item \code{.processSignatureResponseError} for HTTP error responses
#'   \item \code{.processSuccessfulResponse} for successful responses with data
#' }
#'
#' The resulting tibble contains these columns:
#' \itemize{
#'   \item \code{signatureID}: The signature identifier
#'   \item \code{ID_geneid}: Character gene identifiers
#'   \item \code{Name_GeneSymbol}: Gene symbols
#'   \item \code{Value_LogDiffExp}: Log fold-change values (rounded to 12 decimal places)
#'   \item \code{Significance_pvalue}: P-values (rounded to 12 decimal places)
#' }
#'
#' @importFrom httr2 resp_status
#'
#' @keywords internal
.processSignatureResponse <- function(response) {
    # Check response status first
    status <- httr2::resp_status(response)
    if (status != 200L) {
        .processSignatureResponseError(response)
    }

    # Process successful response
    .processSuccessfulResponse(response)
}

#' Get the L1000 Signature from iLINCS
#' `r lifecycle::badge("stable")`
#'
#' This function acts as the entrypoint to the iLINCS database.
#' This takes in an ID and returns the signature after making a
#' call to the iLINCS database. The function automatically detects
#' whether the signature is an L1000 signature based on the signature ID
#' and metadata tables, and retrieves all available genes for comprehensive
#' signature analysis.
#'
#' @param sigId character. The ilincs signature_id
#'
#' @return a tibble with the signature data containing the following columns:
#'   \itemize{
#'     \item \code{signatureID}: The signature identifier
#'     \item \code{ID_geneid}: Gene IDs (Entrez)
#'     \item \code{Name_GeneSymbol}: Gene symbols
#'     \item \code{Value_LogDiffExp}: Log fold-change values
#'     \item \code{Significance_pvalue}: Statistical significance p-values
#'   }
#' @export
#'
#' @importFrom httr2 request req_method req_url_query req_user_agent
#' @importFrom httr2 req_url_path_append req_perform resp_status
#' @importFrom httr2 resp_body_json resp_body_string req_error req_options
#' @importFrom tibble tibble as_tibble
#' @importFrom rlang .data
#' @importFrom dplyr select mutate bind_rows
#' @importFrom purrr map flatten_dfr
#' @importFrom S4Vectors DataFrame
#' @importFrom utils packageVersion
#'
#' @examples
#' # Input validation example (no API call)
#' # Demonstrates proper signature ID format validation
#' tryCatch(
#'     getSignature(""), # Empty string should error
#'     error = function(e) message("Expected error: empty signature ID")
#' )
#'
#' \donttest{
#' # These examples require network access to the iLINCS API
#'
#' # Get the L1000 signature for LINCSKD_28
#' kdSignature <- getSignature("LINCSKD_28")
#' head(kdSignature)
#'
#' # Get an overexpression signature (L1000 status is automatically detected)
#' oeSignature <- getSignature("LINCSOE_1000")
#' head(oeSignature)
#'
#' # Check the structure of retrieved signature
#' str(kdSignature)
#' }
getSignature <- function(sigId) {
    # Validate input parameters
    .validateGetSignatureInput(sigId)

    # Create and execute HTTP request with error handling
    .createSignatureRequest(sigId) |>
        .executeSignatureRequest() |>
        .processSignatureResponse()
}
