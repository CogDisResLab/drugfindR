#' Validate getConcordants input parameters
#'
#' This internal function validates the input parameters for the getConcordants
#' function to ensure they meet the required constraints.
#'
#' @param signature A data.frame-like object containing the signature data.
#'   Must have the required iLINCS signature structure with columns:
#'   signatureID, ID_geneid, Name_GeneSymbol, Value_LogDiffExp, Significance_pvalue.
#' @param ilincsLibrary Character string specifying the iLINCS library to search.
#'   Must be one of "OE", "KD", or "CP".
#'
#' @return A character vector containing the class of the input signature.
#'   This is used internally to determine the return type format.
#'
#' @details
#' This function performs comprehensive validation:
#' \enumerate{
#'   \item Ensures \code{signature} is a data.frame-like object (data.frame, tibble, or S4Vectors::DataFrame)
#'   \item Validates complete signature structure via \code{stopIfInvalidSignature()}
#'   \item Validates \code{ilincsLibrary} is one of the supported libraries
#' }
#'
#' The signature must conform to the iLINCS expected structure. Use
#' \code{prepareSignature()} to ensure proper formatting.
#'
#' @seealso
#' \code{\link{stopIfInvalidSignature}} for signature structure validation,
#' \code{\link{prepareSignature}} for signature preparation
#'
#' @keywords internal
#'
#' @examples NULL
.validateGetConcordantsInput <- function(signature, ilincsLibrary) {
    # 1. Validate signature data type
    if (!any(c("data.frame", "DFrame") %in% class(signature))) {
        stop("Signature must be a data.frame, tibble, or DataFrame", call. = FALSE)
    }

    # 2. Validate the signature structure
    stopIfInvalidSignature(signature)

    # 3. Validate iLINCS library
    stopIfInvalidLibraries(ilincsLibrary)

    invisible(class(signature))
}

#' Prepare signature file for iLINCS upload
#'
#' This internal function creates a temporary file containing the signature data
#' formatted for upload to the iLINCS API.
#'
#' @param signature A data.frame-like object containing the signature data.
#'
#' @return Character string path to the temporary signature file.
#'
#' @details
#' The function creates a temporary file with a ".xls" extension and writes
#' the signature data as a tab-separated file. The file is automatically
#' cleaned up by the system when the R session ends.
#'
#' @importFrom readr write_tsv
#'
#' @keywords internal
#'
#' @examples NULL
.prepareSignatureFile <- function(signature) {
    signatureFile <- tempfile(pattern = "ilincs_sig", fileext = ".xls")
    if ("DFrame" %in% class(signature)) {
        signature <- as.data.frame(signature)
    }
    readr::write_tsv(signature, signatureFile)
    signatureFile
}

#' Detect signature direction from expression values
#'
#' This internal function analyzes the log fold-change values in a signature
#' to determine the overall direction of regulation.
#'
#' @param signature A data.frame-like object containing the signature data.
#'   Must have a column named "Value_LogDiffExp" with log fold-change values.
#'
#' @return Character string indicating signature direction: "Up", "Down", or "Any".
#'
#' @details
#' The function examines the "Value_LogDiffExp" column to determine direction:
#' \itemize{
#'   \item "Up": All expression values are greater than or equal to zero
#'   \item "Down": All expression values are less than or equal to zero
#'   \item "Any": Mixed positive and negative values
#' }
#'
#' Note that zero values are considered "Up" direction. This direction
#' information is used by iLINCS for signature analysis and is included
#' in the output results.
#'
#' @keywords internal
#'
#' @examples NULL
.detectSignatureDirection <- function(signature) {
    expressionValues <- signature[["Value_LogDiffExp"]]

    if (all(expressionValues >= 0L)) {
        "Up"
    } else if (all(expressionValues <= 0L)) {
        "Down"
    } else {
        "Any"
    }
}

#' Create iLINCS API request
#'
#' This internal function constructs and executes the HTTP request to the
#' iLINCS API for concordant signature analysis.
#'
#' @param signatureFile Character string path to the signature file to upload.
#' @param ilincsLibrary Character string specifying the iLINCS library to search.
#'   Must be one of "OE", "KD", or "CP".
#'
#' @return An httr2 response object from the iLINCS API.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Maps the library name to the internal iLINCS library ID
#'   \item Constructs a multipart POST request with the signature file
#'   \item Includes appropriate user agent and API endpoint
#'   \item Executes the request and returns the response
#' }
#'
#' The library mapping is:
#' \itemize{
#'   \item CP (Chemical Perturbagen): LIB_5
#'   \item KD (Knockdown): LIB_6
#'   \item OE (Overexpression): LIB_11
#' }
#'
#' @importFrom httr2 request req_url_path_append req_url_query req_body_multipart
#'   req_method req_user_agent req_error req_perform req_options
#' @importFrom curl form_file
#'
#' @keywords internal
#'
#' @examples NULL
.generateIlincsRequest <- function(signatureFile, ilincsLibrary) {
    # Map library names to iLINCS internal IDs

    # Construct the request
    httr2::request(.ilincsBaseUrl()) |>
        httr2::req_url_path_append("SignatureMeta") |>
        httr2::req_url_path_append("uploadAndAnalyze") |>
        httr2::req_url_query(lib = .returnLibrary(ilincsLibrary)) |>
        httr2::req_body_multipart(file = curl::form_file(signatureFile)) |>
        httr2::req_method("POST") |>
        httr2::req_user_agent(.returnUserAgent())
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
#' processed by \code{.processIlincsResponse()} to generate appropriate
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
#' @importFrom httr2 req_error req_perform
#'
#' @keywords internal
#'
#' @examples NULL
.executeIlincsRequest <- function(request, verbose = FALSE) {
    # Configure request to not raise errors on HTTP error status codes
    request <- httr2::req_error(request, is_error = function(resp) FALSE)

    # Optionally add verbosity for debugging
    if (verbose) {
        request <- httr2::req_options(request, verbose = TRUE)
    }

    # Execute the request - errors are captured in response object
    httr2::req_perform(request)
}

#' Handle iLINCS API response errors
#'
#' This internal function processes error responses from the iLINCS API
#' and generates appropriate error messages.
#'
#' @param response An httr2 response object from the iLINCS API.
#'
#' @return This function always stops execution with an error message.
#'
#' @keywords internal
.processIlincsResponseError <- function(response) {
    status <- httr2::resp_status(response)

    # Extract error details from response body
    errorBody <- tryCatch(
        httr2::resp_body_string(response),
        error = function(e) "Unable to retrieve error details"
    )

    # Create error message based on status code
    if (status == 400L) {
        errorMsg <- paste0(
            "iLINCS API request failed (Status 400): Bad Request\n",
            "This typically indicates invalid signature data or parameters.\n",
            "Ensure the signature has the correct structure using `prepareSignature()`.\n",
            "Response: ", errorBody
        )
    } else if (status == 500L) {
        errorMsg <- paste0(
            "iLINCS API request failed (Status 500): Internal Server Error\n",
            "This may indicate data structure problems or server issues.\n",
            "Verify signature format and try again.\n",
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

#' Handle empty concordance table responses
#'
#' This internal function creates an empty tibble with the correct structure
#' when the iLINCS API returns no concordant signatures.
#'
#' @param sigDirection Character string indicating the signature direction.
#' @param ilincsLibrary Character string specifying the iLINCS library used.
#'
#' @return A tibble with zero rows and the correct column structure.
#'
#' @keywords internal
.processIlincsResponseEmpty <- function(sigDirection, ilincsLibrary) {
    # Define sig_type based on library
    sigType <- switch(ilincsLibrary,
        CP = "Chemical Perturbagen",
        KD = "Gene Knockdown",
        OE = "Gene Overexpression"
    )

    # Return empty tibble with correct column structure
    tibble::tibble(
        signatureid = character(0L),
        treatment = character(0L),
        concentration = character(0L),
        time = character(0L),
        cellline = character(0L),
        similarity = numeric(0L),
        pValue = numeric(0L),
        sig_direction = character(0L),
        sig_type = character(0L)
    )
}

#' Handle successful concordance table responses
#'
#' This internal function processes successful responses from the iLINCS API
#' and formats the concordant signature data with standardized columns.
#'
#' @param concordanceTables List containing concordance table data from API response.
#' @param sigDirection Character string indicating the signature direction.
#' @param ilincsLibrary Character string specifying the iLINCS library used.
#'
#' @return A tibble containing processed concordant signature data.
#'
#' @keywords internal
.processIlincsResponseSuccess <- function(concordanceTables, sigDirection, ilincsLibrary) {
    # Define sig_type based on library
    sigType <- switch(ilincsLibrary,
        CP = "Chemical Perturbagen",
        KD = "Gene Knockdown",
        OE = "Gene Overexpression"
    )

    # Process the concordance tables
    concordants <- concordanceTables |>
        dplyr::bind_rows() |>
        dplyr::select(dplyr::any_of(c(
            "signatureid", "compound", "treatment",
            "concentration", "time", "cellline", "similarity", "pValue"
        )))

    # Handle library-specific column transformations
    if (ilincsLibrary == "CP") {
        # For CP library: rename compound to treatment if it exists
        if ("compound" %in% colnames(concordants)) {
            concordants <- dplyr::rename(concordants, treatment = !!"compound")
        }
    } else {
        # For KD/OE libraries: add concentration column with all NA values
        concordants <- dplyr::mutate(concordants, concentration = NA_character_)
    }

    # Apply final formatting and add metadata columns
    concordants <- concordants |>
        dplyr::select(dplyr::any_of(c(
            "signatureid", "treatment", "concentration",
            "time", "cellline", "similarity", "pValue"
        ))) |>
        dplyr::mutate(
            similarity = round(.data[["similarity"]], 12L),
            pValue = round(.data[["pValue"]], 20L),
            sig_direction = sigDirection,
            sig_type = sigType
        ) |>
        dplyr::select(c(
            "signatureid", "treatment", "concentration", "time",
            "cellline", "similarity", "pValue", "sig_direction", "sig_type"
        ))

    concordants
}

#' Process iLINCS API response into concordant signatures
#'
#' This internal function dispatches to appropriate handlers based on the
#' response status and content from the iLINCS API.
#'
#' @param response An httr2 response object from the iLINCS API.
#' @param sigDirection Character string indicating the signature direction
#'   ("Up", "Down", or "Any").
#' @param ilincsLibrary Character string specifying the iLINCS library used
#'   ("CP", "KD", or "OE").
#'
#' @return A tibble containing concordant signature data with standardized
#'   column names and rounded numerical values.
#'
#' @details
#' The function dispatches to specialized handlers:
#' \enumerate{
#'   \item \code{.processIlincsResponseError} for HTTP error responses
#'   \item \code{.processIlincsResponseEmpty} for empty concordance tables
#'   \item \code{.processIlincsResponseSuccess} for successful responses with data
#' }
#'
#' The resulting tibble always contains these columns in order:
#' \itemize{
#'   \item \code{signatureid}: Unique signature identifier
#'   \item \code{treatment}: Drug/treatment name (compound renamed for CP library)
#'   \item \code{concentration}: Drug concentration (NA for KD/OE libraries)
#'   \item \code{time}: Treatment duration
#'   \item \code{cellline}: Cell line used
#'   \item \code{similarity}: Similarity score (rounded to 8 decimal places)
#'   \item \code{pValue}: Statistical significance (rounded to 20 decimal places)
#'   \item \code{sig_direction}: Signature direction ("Up", "Down", or "Any")
#'   \item \code{sig_type}: Library type description
#' }
#'
#' @importFrom httr2 resp_status resp_body_json
#' @importFrom purrr pluck
#' @importFrom dplyr bind_rows select any_of rename mutate
#' @importFrom rlang .data
#'
#' @keywords internal
.processIlincsResponse <- function(response, sigDirection, ilincsLibrary) {
    # Check response status first
    status <- httr2::resp_status(response)
    if (status != 200L) {
        .processIlincsResponseError(response)
    }

    # Process successful response
    responseData <- httr2::resp_body_json(response)
    concordanceTables <- purrr::pluck(responseData, "status", "concordanceTable")

    # Check if concordanceTable is empty
    # The API returns an empty list in case of error
    if (length(concordanceTables) == 0L) {
        .processIlincsResponseEmpty(sigDirection, ilincsLibrary)
    } else {
        .processIlincsResponseSuccess(concordanceTables, sigDirection, ilincsLibrary)
    }
}


#' Clean up temporary signature file
#'
#' This internal function removes the temporary signature file created during
#' the getConcordants operation to prevent accumulation of temporary files.
#'
#' @param signatureFile Character string path to the temporary signature file
#'   to be removed.
#'
#' @return Invisible NULL. The function is called for its side effect of
#'   removing the temporary file.
#'
#' @details
#' The function checks if the specified file exists and removes it using
#' \code{unlink()}. This cleanup is performed automatically at the end of
#' the getConcordants operation.
#'
#' @keywords internal
#'
#' @examples NULL
.cleanupGetConcordants <- function(signatureFile) {
    if (file.exists(signatureFile)) unlink(signatureFile)
}


#' Return results in appropriate format based on input type
#'
#' This internal function formats the output results to match the input
#' signature type, ensuring consistent data type handling.
#'
#' @param result A tibble containing the processed results from iLINCS API.
#' @param inputClass A character vector containing the class of the original
#'   input signature (from \code{.validateGetConcordantsInput}).
#'
#' @return The results in the appropriate format:
#'   \itemize{
#'     \item S4Vectors::DataFrame if input was a DataFrame
#'     \item tibble otherwise (for data.frame, tibble inputs)
#'   }
#'
#' @details
#' This function ensures that the output format matches the input format
#' for consistency. If the original signature was provided as an
#' S4Vectors::DataFrame, the results are converted back to DataFrame.
#' Otherwise, results are returned as a tibble.
#'
#' @importFrom S4Vectors DataFrame
#'
#' @keywords internal
#'
#' @examples NULL
.returnResults <- function(result, inputClass) {
    if ("DFrame" %in% inputClass) {
        DataFrame(result)
    } else {
        result
    }
}

#' Get concordant signatures from iLINCS database
#'
#' This function queries the iLINCS (Integrative Library of Integrated
#' Network-based Cellular Signatures) database to find signatures that are
#' concordant (similar) to a given input signature.
#'
#' @param signature A data.frame, tibble, or S4Vectors::DataFrame containing
#'   the signature data. Must conform to iLINCS signature structure with columns:
#'   \itemize{
#'     \item \code{signatureID}: Signature identifier
#'     \item \code{ID_geneid}: Gene IDs
#'     \item \code{Name_GeneSymbol}: Gene symbols
#'     \item \code{Value_LogDiffExp}: Log fold-change values
#'     \item \code{Significance_pvalue}: Statistical significance p-values
#'   }
#'   Use \code{prepareSignature()} to ensure proper formatting.
#' @param ilincsLibrary Character string specifying the iLINCS library to search.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"CP"}: Chemical Perturbagen library (default)
#'     \item \code{"KD"}: Knockdown library
#'     \item \code{"OE"}: Overexpression library
#'   }
#'
#' @return A data structure containing concordant signatures. The return type
#'   matches the input signature type:
#'   \itemize{
#'     \item tibble for data.frame or tibble inputs
#'     \item S4Vectors::DataFrame for DataFrame inputs
#'   }
#'
#'   Contains the following columns:
#'   \itemize{
#'     \item \code{signatureid}: Unique signature identifier
#'     \item \code{compound} or \code{treatment}: Drug/treatment name
#'     \item \code{concentration}: Drug concentration (CP library only)
#'     \item \code{time}: Treatment duration
#'     \item \code{cellline}: Cell line used
#'     \item \code{similarity}: Similarity score (rounded to 8 decimal places)
#'     \item \code{pValue}: Statistical significance (rounded to 20 decimal places)
#'     \item \code{sig_direction}: Signature direction ("Up", "Down", or "Any")
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Creates a temporary file with signature data
#'   \item Detects signature direction from expression values
#'   \item Sends a multipart POST request to the iLINCS API
#'   \item Processes the JSON response into a standardized tibble
#'   \item Cleans up temporary files
#' }
#'
#' The signature direction is determined as follows:
#' \itemize{
#'   \item \code{"Up"}: All expression values are greater than or equal to zero
#'   \item \code{"Down"}: All expression values are less than or equal to zero
#'   \item \code{"Any"}: Mixed positive and negative values
#' }
#'
#' @section API Details:
#' This function interfaces with the iLINCS web service API. The signature
#' is uploaded as a tab-separated file and analyzed against the specified
#' library. Results are returned as JSON and parsed into a tibble.
#'
#' @section Error Handling:
#' The function will stop execution with informative error messages for:
#' \itemize{
#'   \item Invalid signature data types (must be data.frame, tibble, or DataFrame)
#'   \item Invalid signature structure (missing required columns, wrong order, etc.)
#'   \item Missing values in signature data
#'   \item Unsupported iLINCS library names
#'   \item HTTP errors from the iLINCS API
#'   \item Invalid or empty API responses
#' }
#'
#' @importFrom httr2 request req_url_path_append req_url_query req_body_multipart
#'   req_method req_user_agent req_error req_perform req_options resp_status resp_body_json resp_body_string
#' @importFrom curl form_file
#' @importFrom readr write_tsv
#' @importFrom purrr map flatten_dfr
#' @importFrom dplyr select any_of mutate
#' @importFrom rlang .data
#' @importFrom S4Vectors DataFrame
#'
#' @seealso
#' \code{\link{prepareSignature}} for signature preparation,
#' \code{\link{filterSignature}} for signature filtering,
#' \code{\link{investigateSignature}} for signature investigation
#'
#' @references
#' iLINCS Portal: \url{http://www.ilincs.org/}
#'
#' Pilarczyk et al. (2020). Connecting omics signatures and revealing
#' biological mechanisms with iLINCS. Nature Communications, 11(1), 4058.
#'
#' @examples
#' # Input validation examples (no API calls)
#' # These demonstrate proper signature structure
#' mockSig <- data.frame(
#'     signatureID = rep("TEST", 3),
#'     ID_geneid = c("123", "456", "789"),
#'     Name_GeneSymbol = c("TP53", "MYC", "EGFR"),
#'     Value_LogDiffExp = c(1.5, -2.0, 0.8)
#' )
#'
#' # Validate library parameter (should produce error)
#' tryCatch(
#'     getConcordants(mockSig, ilincsLibrary = "INVALID"),
#'     error = function(e) message("Expected error: invalid library")
#' )
#'
#' \donttest{
#' # This example requires network access to the iLINCS API
#'
#' # Load example differential expression data
#' dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
#'     package = "drugfindR"
#' )
#' dge_data <- read.delim(dge_file)
#'
#' # Prepare signature to ensure proper structure
#' signature <- prepareSignature(
#'     dge_data[1:50, ],
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC",
#'     pvalColumn = "PValue"
#' )
#'
#' # Find concordant chemical perturbagens
#' cpConcordants <- getConcordants(signature, ilincsLibrary = "CP")
#' head(cpConcordants)
#'
#' # Find concordant knockdown signatures
#' kdConcordants <- getConcordants(signature, ilincsLibrary = "KD")
#' head(kdConcordants)
#'
#' # Find concordant overexpression signatures
#' oeConcordants <- getConcordants(signature, ilincsLibrary = "OE")
#' head(oeConcordants)
#'
#' # Works with different data frame types
#' signatureDf <- as.data.frame(signature)
#' cpConcordantsDf <- getConcordants(signatureDf, "CP")
#'
#' # Works with S4Vectors::DataFrame
#' signatureDataFrame <- S4Vectors::DataFrame(signature)
#' cpConcordantsDataFrame <- getConcordants(signatureDataFrame, "CP")
#' # Returns S4Vectors::DataFrame to match input type
#' }
#'
#' @export
getConcordants <- function(signature, ilincsLibrary = "CP") {
    # Validate input parameters
    inputClass <- .validateGetConcordantsInput(signature, ilincsLibrary)

    # Prepare signature file for upload
    signatureFile <- .prepareSignatureFile(signature)

    # Detect signature direction
    sigDirection <- .detectSignatureDirection(signature)

    # Generate iLINCS API request
    request <- .generateIlincsRequest(signatureFile, ilincsLibrary)

    # Execute the request with error handling
    response <- .executeIlincsRequest(request)

    # Process response and return results
    result <- .processIlincsResponse(response, sigDirection, ilincsLibrary)

    # Cleanup the generated files
    .cleanupGetConcordants(signatureFile)

    # Return the results
    .returnResults(result, inputClass)
}
