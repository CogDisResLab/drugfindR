#' @include class-drugfindRDataset.R
NULL

#' Construct a drugfindRDataset Object
#'
#' This function wraps a `tibble`, `data.frame`, or `DataFrame` into a
#' `drugfindRDataset` S4 object, for standardized manipulation in the
#' `drugfindR` pipeline.
#'
#' @param dflike A data-like object — `data.frame`, `tibble`, or `DataFrame`
#'
#' @return An S4 object of class `drugfindRDataset`
#' @export
#'
#' @importFrom methods new
#'
#' @examples
#' df <- tibble::tibble(Gene = letters[1:5], Value_LogDiffExp = c(-1, -0.5, 0, 0.5, 1))
#' sig <- drugfindRDataset(df)
#' sig
drugfindRDataset <- function(dflike) {
    if (!inherits(dflike, c("data.frame", "tbl_df", "DFrame"))) {
        stop("Input must be a data.frame, tibble, or S4Vectors::DataFrame.")
    }

    df <- as_tibble(dflike)
    core <- drugfindRCoreData(signature = df)
    source_class <- class(dflike)[1]

    methods::new("drugfindRDataset", core = core, source_class = source_class)
}
