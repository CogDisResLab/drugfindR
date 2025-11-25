#' Rename target-related columns to user-facing output names
#'
#' Standardizes internal concordants/consensus column names to the
#' user-facing output schema used by this package.
#'
#' @param inputNames Character vector of column names to rename. See Details
#'   for the expected input ordering and mapping.
#'
#' @return Character vector of output (renamed) column names.
#'
#' @details
#' Expected input columns (by position) are the internal concordants fields:
#'
#' 1) `signatureid`, 2) `treatment` (or `compound` pre-renaming),
#' 3) `cellline`, 4) `time`, 5) `concentration`, 6) `sig_direction`,
#' 7) `sig_type`, 8) `similarity`, 9) `pValue`.
#'
#' These are mapped to the user-facing names returned by functions like
#' [consensusConcordants()] and downstream investigation helpers:
#'
#' - `TargetSignature`, `Target`, `TargetCellLine`, `TargetTime`,
#'   `TargetConcentration`, `InputSigDirection`, `SignatureType`,
#'   `Similarity`, `pValue`.
#'
#' Only the names are returned; the renaming is applied via
#' `dplyr::rename_with(x, targetRename)`.
#'
#' @keywords internal
#'
#' @examples NULL
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
#' @examples NULL
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
#' @examples NULL
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
#' @examples NULL
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
#' @examples NULL
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
            errorMsg <- paste0(errorMsg, "Missing columns: ", toString(missingCols), "\n")
        }

        if (length(extraCols) > 0L) {
            errorMsg <- paste0(errorMsg, "Unexpected columns: ", toString(extraCols), "\n")
        }

        if (!identical(expectedColNames, actualColNames[seq_along(expectedColNames)])) {
            errorMsg <- paste0(errorMsg, "Columns are not in the expected order.\n")
        }

        errorMsg <- paste0(
            errorMsg,
            "Expected columns (in order): ", toString(expectedColNames), "\n",
            "Actual columns: ", toString(actualColNames), "\n",
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
#' @examples NULL
.stopIfContainsMissingValues <- function(signature) {
    if (anyNA(signature)) {
        # Find which columns contain missing values
        colsWithNa <- colnames(signature)[vapply(signature, anyNA, logical(1L))]
        naCounts <- vapply(signature[colsWithNa], function(x) sum(is.na(x)), integer(1L))

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
#' @examples NULL
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
#' @examples NULL
.loadMetadata <- function(lib) {
    switch(lib,
        OE = oeMetadata,
        KD = kdMetadata,
        CP = cpMetadata,
        stop(
            "Invalid library: '", lib, "'. ",
            "Library must be one of 'OE' (Overexpression), 'KD' (Knockdown), or 'CP' (Chemical Perturbagen).",
            call. = FALSE
        )
    )
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
#' @examples NULL
.returnLibrary <- function(lib) {
    stopIfInvalidLibraries(lib)

    switch(lib,
        OE = "LIB_11",
        KD = "LIB_6",
        CP = "LIB_5"
    )
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
#' @examples NULL
.returnUserAgent <- function() {
    paste0(
        "drugfindR/",
        utils::packageVersion("drugfindR"),
        "; https://github.com/CogDisResLab/drugfindR"
    )
}

#' Compute consensus concordant signatures from a single input signature
#'
#' This internal helper wraps the common paired / unpaired workflow used by
#' `investigateSignature()` and `investigateTarget()` for a single already
#' prepared or retrieved signature. It applies directional filtering, queries
#' iLINCS for concordant signatures, and collapses results via
#' [consensusConcordants()].
#'
#' @param signature A data.frame / tibble / DataFrame produced by
#'   [prepareSignature()] or [getSignature()] with standard signature columns.
#' @param outputLib Character. One of "OE", "KD", or "CP" indicating the
#'   iLINCS library to search for concordant signatures.
#' @param filterThreshold Numeric (optional). Absolute threshold(s) passed to
#'   [filterSignature()]. Use either `filterThreshold` or `filterProp`.
#' @param filterProp Numeric in (0, 0.5] (optional). Proportion for quantile
#'   based filtering in [filterSignature()]. Ignored if `filterThreshold` is
#'   supplied.
#' @param similarityThreshold Numeric in 0..1. Minimum absolute similarity
#'   retained by [consensusConcordants()].
#' @param paired Logical. If TRUE perform separate up / down filtering and
#'   concordance queries; otherwise aggregate direction = "any".
#' @param outputCellLines Optional character vector restricting target cell
#'   lines during consensus filtering. Passed to [consensusConcordants()].
#'
#' @return A tibble of consensus concordant signatures with standardized target
#'   columns (already renamed via internal consensus pipeline). Columns include
#'   `TargetSignature`, `Target`, `TargetCellLine`, `Similarity`, `pValue`,
#'   `InputSigDirection`, `SignatureType`, and optional time / concentration.
#'
#' @details
#' Error handling is delegated to component functions:
#' * Library validation via [stopIfInvalidLibraries()]
#' * Signature structure via [stopIfInvalidSignature()] (indirectly used by
#'   [getConcordants()])
#' * Filtering parameter validation via [.validateFilterSignatureInput()]
#' * Concordance / network errors via internal iLINCS response processors.
#'
#' If both `filterThreshold` and `filterProp` are supplied an error is raised
#' upstream in [filterSignature()]. Provide only one.
#'
#' @keywords internal
#'
#' @importFrom dplyr mutate
#'
#' @examples NULL
.computeConsensusFromSignature <- function(
  signature,
  outputLib,
  filterThreshold = NULL,
  filterProp = NULL,
  similarityThreshold = 0.321,
  paired = TRUE,
  outputCellLines = NULL
) {
    stopIfInvalidLibraries(outputLib)

    if (paired) {
        upSig <- filterSignature(signature,
            direction = "up",
            threshold = filterThreshold, prop = filterProp
        )
        downSig <- filterSignature(signature,
            direction = "down",
            threshold = filterThreshold, prop = filterProp
        )
        upConc <- getConcordants(upSig, ilincsLibrary = outputLib)
        downConc <- getConcordants(downSig, ilincsLibrary = outputLib)
        consensusConcordants(upConc, downConc,
            paired = TRUE,
            cellLine = outputCellLines,
            cutoff = similarityThreshold
        )
    } else {
        anySig <- filterSignature(signature,
            direction = "any",
            threshold = filterThreshold, prop = filterProp
        )
        conc <- getConcordants(anySig, ilincsLibrary = outputLib)
        consensusConcordants(conc,
            paired = FALSE,
            cellLine = outputCellLines,
            cutoff = similarityThreshold
        )
    }
}
