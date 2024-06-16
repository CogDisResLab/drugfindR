#' Rename the Target-Related Columns
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
