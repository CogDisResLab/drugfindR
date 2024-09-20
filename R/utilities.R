#' Rename the Target-Related Columns
#'
#' This function is used to standardize the names
#' of the columns output at the end of the result.
#'
#' @param inputNames A character vector of input_names
#'
#' @return A character vector of new names
targetRename <- function(inputNames) {
    if ("treatment" %in% inputNames) {
        c(
            "TargetSignature", "Target", "TargetCellLine",
            "TargetTime", "Similarity", "SignatureDirection", "pValue"
        )
    } else {
        c(
            "TargetSignature", "Target", "TargetCellLine",
            "TargetTime", "TargetConcentration", "Similarity",
            "SignatureDirection", "pValue"
        )
    }
}

#' Parameterize the base URL for the iLINCS API
#'
#'
#' @return a fixed string URL
.ilincsBaseUrl <- function() {
    "http://www.ilincs.org/api"
}

#' Check if the library is valid
#'
#' @param lib a string of libraries
#'
#' @return a boolean
.validateLibrary <- function(lib) {
    lib %in% c("CP", "KD", "OE")
}

#' Check if the libraries input are valid
#'
#' This function confirms that the library parameter
#' is one of the allowed ones.
#'
#' @param libs a character vector of libraries
#'
#' @return a boolean
validateLibraries <- function(libs) {
    all(purrr::map_lgl(libs, .validateLibrary))
}

#' Stop if the libraries are invalid
#'
#' @param libs a character vector of libraries
#'
#' @return a stop if the libraries are invalid
stopIfInvalidLibraries <- function(libs) {
    if (!validateLibraries(libs)) {
        stop("Both input and output libraries must be one of 'OE', 'KD', 'CP'")
    }
}

#' Load the correct metadata table
#'
#' @param lib a string. One of "OE", "KD" or "CP"
#'
#' @return a tibble
loadMetadata <- function(lib) {
    if (lib == "OE") {
        oeMetadata
    } else if (lib == "KD") {
        kdMetadata
    } else if (lib == "CP") {
        cpMetadata
    } else {
        stop("Invalid library")
    }
}


#' Return the internal iLINCS Library ID for a given library
#'
#' @param lib A library name. Can be one of "OE", "KD" or "CP"
#'
#' @return A string with the associated library ID
.return_library <- function(lib) {
    stopIfInvalidLibraries(lib)

    libMap <- c(
        OE = "LIB_11",
        KD = "LIB_6",
        CP = "LIB_5"
    )

    libMap[lib]
}
