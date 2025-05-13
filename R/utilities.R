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
#' @keywords internal
#' @return a fixed string URL
.ilincsBaseUrl <- function() {
    "https://www.ilincs.org/api"
}

#' Check if the library is valid
#'
#' @param lib a string of libraries
#'
#' @keywords internal
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
#' @keywords internal
#' @return a boolean
validateLibraries <- function(libs) {
    all(purrr::map_lgl(libs, .validateLibrary))
}

#' Stop if the libraries are invalid
#'
#' @param libs a character vector of libraries
#'
#' @keywords internal
#' @return a stop if the libraries are invalid
stopIfInvalidLibraries <- function(libs) {
    if (!validateLibraries(libs)) {
        stop("Both input and output libraries must be one of 'OE', 'KD', 'CP'")
    }
}

#' Check if the cell line is valid
#'
#' @param cellLine a string of cell lines
#'
#' @keywords internal
#' @return a boolean
.validateCellLine <- function(cellLine) {
    validCellLines <- c(
        "A375", "HEK293T", "HEPG2", "HT29", "MCF7", "PC3", "VCAP",
        "A549", "HA1E", "HCC515", "HELA", "A673", "AGS", "ASC", "HME1",
        "HS578T", "BT20", "CD34", "CL34", "CORL23", "COV644", "DV90",
        "HUES3", "HUVEC", "EFO27", "JURKAT", "LNCAP", "FIBRNPC", "MCF10A",
        "H1299", "MDAMB231", "MNEU.E", "NEU", "NPC.CAS9", "NPC", "NPC.TAK",
        "SKBR3", "SKL.C", "SKL", "YAPC", "HCC15", "HCT116", "HEC108",
        "HL60", "HT115", "HUH7", "JHUEM2", "LOVO", "HS27A", "MDST8",
        "NCIH1694", "NCIH1836", "NCIH2073", "NCIH508", "NCIH596", "NKDBA",
        "NOMO1", "OV7", "PHH", "PL21", "RKO", "RMGI", "RMUGS", "SKB",
        "SKLU1", "SKM1", "SKMEL1", "SKMEL28", "SNGM", "SNU1040", "SNUC4",
        "SNUC5", "SW480", "SW620", "SW948", "T3M10", "THP1", "TYKNU",
        "U266", "U937", "WSUDLCL2", "MCH58", "ASC.C", "HEKTE", "NCIH716",
        "SHSY5Y", "A375.311", "A549.311", "HA1E.101", "HA1E.311", "HELA.311",
        "HT29.311", "MCF7.101", "MCF7.311", "PC3.101", "PC3.311", "YAPC.311"
    )
    cellLine %in% validCellLines
}

#' Check if the cell lines input are valid
#'
#' This function confirms that the library parameter
#' is one of the allowed ones.
#'
#' @param cellLines a character vector of libraries
#'
#' @keywords internal
#' @return a boolean
validateCellLines <- function(cellLines) {
    all(purrr::map_lgl(cellLines, .validateCellLine))
}

#' Load the correct metadata table
#'
#' @param lib a string. One of "OE", "KD" or "CP"
#'
#' @keywords internal
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
#' @keywords internal
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

#' Return a string suitable as a User-Agent for the iLINCS API
#'
#' @keywords internal
#' @return a string
.return_user_agent <- function() {
    paste0("drugfindR/", utils::packageVersion("drugfindR"), "; https://github.com/CogDisResLab/drugfindR")
}
