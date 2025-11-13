#' Rename the Target-Related Columns
#'
#' This function is used to standardize the names
#' of the columns output at the end of the result.
#'
#' @param inputNames A character vector of input_names
#'
#' @return A character vector of new names
targetRename <- function(inputNames) {
    c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )
}

#' Parameterize the base URL for the iLINCS API
#'
#' @keywords internal
#' @return a fixed string URL
.ilincsBaseUrl <- function() {
    "https://www.ilincs.org/api"
}

#' Check if a single library is valid
#'
#' This internal function validates whether a single library name is one of
#' the supported iLINCS library types.
#'
#' @param lib A character string containing a single library name to validate.
#'
#' @return A logical value: TRUE if the library is valid, FALSE otherwise.
#'
#' @details
#' Valid library names are:
#' \itemize{
#'   \item \code{"CP"}: Chemical Perturbagen library
#'   \item \code{"KD"}: Knockdown library
#'   \item \code{"OE"}: Overexpression library
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' .validateLibrary("CP") # Returns TRUE
#' .validateLibrary("KD") # Returns TRUE
#' .validateLibrary("OE") # Returns TRUE
#' .validateLibrary("INVALID") # Returns FALSE
#' }
.validateLibrary <- function(lib) {
    lib %in% c("CP", "KD", "OE")
}

#' Check if multiple libraries are valid
#'
#' This function validates whether all provided library names are supported
#' iLINCS library types.
#'
#' @param libs A character vector of library names to validate.
#'
#' @return A logical value: TRUE if all libraries are valid, FALSE if any are invalid.
#'
#' @details
#' This function uses \code{.validateLibrary()} to check each library individually
#' and returns TRUE only if all libraries are valid. It's used internally to
#' validate library parameters before API calls.
#'
#' @importFrom purrr map_lgl
#'
#' @seealso
#' \code{\link{.validateLibrary}} for single library validation,
#' \code{\link{stopIfInvalidLibraries}} for validation with error handling
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' validateLibraries(c("CP", "KD")) # Returns TRUE
#' validateLibraries("OE") # Returns TRUE
#' validateLibraries(c("CP", "INVALID")) # Returns FALSE
#' validateLibraries("INVALID") # Returns FALSE
#' }
validateLibraries <- function(libs) {
    all(purrr::map_lgl(libs, .validateLibrary))
}

#' Stop if the libraries are invalid
#'
#' This internal function validates library specifications and stops execution
#' with an informative error message if any invalid libraries are found.
#'
#' @param libs A character vector of library names to validate.
#'   Each library must be one of "OE", "KD", or "CP".
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' This function validates that all provided library names are supported:
#' \itemize{
#'   \item \code{"OE"}: Overexpression library (LIB_11)
#'   \item \code{"KD"}: Knockdown library (LIB_6)
#'   \item \code{"CP"}: Chemical Perturbagen library (LIB_5)
#' }
#'
#' If any invalid libraries are found, the function provides a detailed error
#' message listing the invalid libraries and the expected options.
#'
#' @importFrom purrr map_lgl
#'
#' @seealso
#' \code{\link{validateLibraries}} for the underlying validation logic,
#' \code{\link{.validateLibrary}} for single library validation
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid libraries
#' stopIfInvalidLibraries(c("CP", "KD")) # No error
#' stopIfInvalidLibraries("OE") # No error
#'
#' # Invalid libraries (will throw errors)
#' stopIfInvalidLibraries("INVALID") # Error
#' stopIfInvalidLibraries(c("CP", "XYZ")) # Error
#' }
stopIfInvalidLibraries <- function(libs) {
    if (!validateLibraries(libs)) {
        invalidLibs <- libs[!purrr::map_lgl(libs, .validateLibrary)]
        stop(
            "Invalid library specification(s): ", paste(invalidLibs, collapse = ", "), ". ",
            "Libraries must be one of 'OE' (Overexpression), 'KD' (Knockdown), or 'CP' (Chemical Perturbagen).",
            call. = FALSE
        )
    }
}

#' Validate signature column names
#'
#' This internal function checks if the signature data frame has the expected
#' column names in the correct order for iLINCS compatibility.
#'
#' @param signature A data.frame-like object containing signature data.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' The function validates that the signature has exactly the following columns
#' in the specified order:
#' \enumerate{
#'   \item \code{signatureID}: Signature identifier
#'   \item \code{ID_geneid}: Gene ID
#'   \item \code{Name_GeneSymbol}: Gene symbol
#'   \item \code{Value_LogDiffExp}: Log fold-change expression value
#'   \item \code{Significance_pvalue}: Statistical significance p-value
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid signature structure
#' valid_sig <- data.frame(
#'     signatureID = "SIG_001",
#'     ID_geneid = "1234",
#'     Name_GeneSymbol = "GENE1",
#'     Value_LogDiffExp = 1.5,
#'     Significance_pvalue = 0.05
#' )
#' .stopIfInvalidColNames(valid_sig) # No error
#'
#' # Invalid signature structure (will throw error)
#' invalid_sig <- data.frame(Gene = "GENE1", Expression = 1.5)
#' .stopIfInvalidColNames(invalid_sig) # Error
#' }
.stopIfInvalidColNames <- function(signature) {
    expectedColNames <- c(
        "signatureID",
        "ID_geneid",
        "Name_GeneSymbol",
        "Value_LogDiffExp",
        "Significance_pvalue"
    )

    actualColNames <- colnames(signature)

    if (!identical(expectedColNames, actualColNames)) {
        missingCols <- setdiff(expectedColNames, actualColNames)
        extraCols <- setdiff(actualColNames, expectedColNames)

        errorMsg <- "Input signature does not conform to expected structure.\n"

        if (length(missingCols) > 0L) {
            errorMsg <- paste0(errorMsg, "Missing columns: ", paste(missingCols, collapse = ", "), "\n")
        }

        if (length(extraCols) > 0L) {
            errorMsg <- paste0(errorMsg, "Unexpected columns: ", paste(extraCols, collapse = ", "), "\n")
        }

        if (!identical(expectedColNames, actualColNames[seq_along(expectedColNames)])) {
            errorMsg <- paste0(errorMsg, "Columns are not in the expected order.\n")
        }

        errorMsg <- paste0(
            errorMsg,
            "Expected columns (in order): ", paste(expectedColNames, collapse = ", "), "\n",
            "Actual columns: ", paste(actualColNames, collapse = ", "), "\n",
            "Please use `prepareSignature()` to ensure compliance."
        )

        stop(errorMsg, call. = FALSE)
    }
}

#' Validate signature for missing values
#'
#' This internal function checks if the signature data frame contains any
#' missing (NA) values, which are not allowed in iLINCS signature data.
#'
#' @param signature A data.frame-like object containing signature data.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' The function scans the entire signature data frame for missing values.
#' iLINCS requires complete data for all signature analysis, so any NA values
#' will cause the function to stop with an informative error message indicating
#' which columns contain missing values.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid signature (no missing values)
#' valid_sig <- data.frame(
#'     signatureID = "SIG_001",
#'     ID_geneid = "1234",
#'     Name_GeneSymbol = "GENE1",
#'     Value_LogDiffExp = 1.5,
#'     Significance_pvalue = 0.05
#' )
#' .stopIfContainsMissingValues(valid_sig) # No error
#'
#' # Invalid signature (contains NA values)
#' invalid_sig <- data.frame(
#'     signatureID = "SIG_001",
#'     ID_geneid = NA,
#'     Name_GeneSymbol = "GENE1",
#'     Value_LogDiffExp = 1.5,
#'     Significance_pvalue = NA
#' )
#' .stopIfContainsMissingValues(invalid_sig) # Error
#' }
.stopIfContainsMissingValues <- function(signature) {
    if (any(is.na(signature))) {
        # Find which columns contain missing values
        colsWithNa <- colnames(signature)[sapply(signature, function(x) any(is.na(x)))]
        naCounts <- sapply(signature[colsWithNa], function(x) sum(is.na(x)))

        errorMsg <- paste0(
            "Input signature contains missing (NA) values, which are not allowed.\n",
            "Columns with missing values:\n"
        )

        for (i in seq_along(colsWithNa)) {
            errorMsg <- paste0(
                errorMsg, "  - ", colsWithNa[i], ": ", naCounts[i], " missing value(s)\n"
            )
        }

        errorMsg <- paste0(
            errorMsg,
            "Please remove or impute missing values before proceeding. ",
            "Consider using `prepareSignature()` with appropriate data cleaning."
        )

        stop(errorMsg, call. = FALSE)
    }
}


#' Validate signature data structure and content
#'
#' This function performs comprehensive validation of signature data to ensure
#' it meets the requirements for iLINCS analysis.
#'
#' @param signature A data.frame-like object containing signature data that
#'   should be validated for iLINCS compatibility.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' This function performs two main validation checks:
#' \enumerate{
#'   \item Column structure validation via \code{.stopIfInvalidColNames()}
#'   \item Missing value validation via \code{.stopIfContainsMissingValues()}
#' }
#'
#' The signature must have exactly the required columns in the correct order
#' and cannot contain any missing (NA) values.
#'
#' @seealso
#' \code{\link{prepareSignature}} for preparing signatures that meet these requirements,
#' \code{\link{.stopIfInvalidColNames}} for column validation details,
#' \code{\link{.stopIfContainsMissingValues}} for missing value validation details
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid signature
#' valid_sig <- data.frame(
#'     signatureID = rep("SIG_001", 3),
#'     ID_geneid = c("1234", "5678", "9012"),
#'     Name_GeneSymbol = c("GENE1", "GENE2", "GENE3"),
#'     Value_LogDiffExp = c(1.5, -2.1, 0.8),
#'     Significance_pvalue = c(0.05, 0.01, 0.03)
#' )
#' stopIfInvalidSignature(valid_sig) # No error
#'
#' # Invalid signature (wrong columns)
#' invalid_sig <- data.frame(Gene = "GENE1", Expression = 1.5)
#' stopIfInvalidSignature(invalid_sig) # Error
#' }
stopIfInvalidSignature <- function(signature) {
    # Ensure that all the required column names are present
    .stopIfInvalidColNames(signature)

    # Ensure that there are no missing values in the signature
    .stopIfContainsMissingValues(signature)
}

#' Load the correct metadata table for a given library
#'
#' This internal function retrieves the appropriate metadata table based on
#' the specified iLINCS library type.
#'
#' @param lib A character string specifying the library type.
#'   Must be one of "OE", "KD", or "CP".
#'
#' @return A tibble containing the metadata for the specified library.
#'   The structure varies by library type but typically includes columns
#'   for signature identifiers, treatments, cell lines, and other metadata.
#'
#' @details
#' The function loads pre-compiled metadata tables for each library:
#' \itemize{
#'   \item \code{"OE"}: Overexpression metadata (\code{oeMetadata})
#'   \item \code{"KD"}: Knockdown metadata (\code{kdMetadata})
#'   \item \code{"CP"}: Chemical Perturbagen metadata (\code{cpMetadata})
#' }
#'
#' These metadata tables are included with the package and contain information
#' about available signatures in each iLINCS library.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Load metadata for different libraries
#' cp_meta <- .loadMetadata("CP")
#' kd_meta <- .loadMetadata("KD")
#' oe_meta <- .loadMetadata("OE")
#'
#' # Invalid library (will throw error)
#' .loadMetadata("INVALID") # Error
#' }
.loadMetadata <- function(lib) {
    if (lib == "OE") {
        oeMetadata
    } else if (lib == "KD") {
        kdMetadata
    } else if (lib == "CP") {
        cpMetadata
    } else {
        stop(
            "Invalid library: '", lib, "'. ",
            "Library must be one of 'OE' (Overexpression), 'KD' (Knockdown), or 'CP' (Chemical Perturbagen).",
            call. = FALSE
        )
    }
}


#' Return the internal iLINCS Library ID for a given library
#'
#' This internal function maps user-friendly library names to the internal
#' library identifiers used by the iLINCS API.
#'
#' @param lib A character string specifying the library name.
#'   Must be one of "OE", "KD", or "CP".
#'
#' @return A character string containing the corresponding iLINCS library ID.
#'
#' @details
#' The mapping between user library names and iLINCS internal IDs is:
#' \itemize{
#'   \item \code{"OE"} (Overexpression) -> \code{"LIB_11"}
#'   \item \code{"KD"} (Knockdown) -> \code{"LIB_6"}
#'   \item \code{"CP"} (Chemical Perturbagen) -> \code{"LIB_5"}
#' }
#'
#' The function validates the input library name before mapping and will
#' stop execution if an invalid library is provided.
#'
#' @seealso
#' \code{\link{stopIfInvalidLibraries}} for library validation details
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid library mappings
#' .returnLibrary("CP") # Returns "LIB_5"
#' .returnLibrary("KD") # Returns "LIB_6"
#' .returnLibrary("OE") # Returns "LIB_11"
#'
#' # Invalid library (will throw error)
#' .returnLibrary("INVALID") # Error via stopIfInvalidLibraries()
#' }
.returnLibrary <- function(lib) {
    stopIfInvalidLibraries(lib)

    libMap <- c(
        OE = "LIB_11",
        KD = "LIB_6",
        CP = "LIB_5"
    )

    libMap[[lib]]
}

#' Return a string suitable as a User-Agent for the iLINCS API
#'
#' This internal function constructs a standardized User-Agent string for
#' HTTP requests to the iLINCS API, including package name, version, and
#' repository URL for identification and debugging purposes.
#'
#' @return A character string formatted as a User-Agent header value.
#'
#' @details
#' The User-Agent string follows the format:
#' \code{"drugfindR/<current version>; https://github.com/CogDisResLab/drugfindR"}
#'
#' This helps iLINCS administrators identify requests from this package
#' and assists with debugging if issues arise. The version is automatically
#' retrieved from the package metadata.
#'
#' @importFrom utils packageVersion
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get the current User-Agent string
#' user_agent <- .returnUserAgent()
#' # Returns something like: "drugfindR/1.0.0; https://github.com/CogDisResLab/drugfindR"
#' }
.returnUserAgent <- function() {
    paste0(
        "drugfindR/",
        utils::packageVersion("drugfindR"),
        "; https://github.com/CogDisResLab/drugfindR"
    )
}
