#' @include class-drugfindRDataset.R
NULL

#' Construct a drugfindRDataset Object
#'
#' This function wraps a `tibble`, `data.frame`, or `DataFrame` into a
#' `drugfindRDataset` S4 object, for standardized manipulation in the
#' `drugfindR` pipeline.
#'
#' @param input A data-like object — `data.frame`, `tibble`, or `DataFrame`
#'
#' @return An S4 object of class `drugfindRDataset`
#' @export
#'
#' @importFrom methods new
#'
#' @examples
#' df <- tibble::tibble(
#'     Gene = letters[1:5],
#'     Value_LogDiffExp = c(-1, -0.5, 0, 0.5, 1)
#' )
#' sig <- drugfindRDataset(df)
#' sig
drugfindRDataset <- function(input) {
    df <- as_tibble(input)
    core <- drugfindRCoreData(coreSignature = df)
    sourceClass <- class(input)[1L]

    methods::new("drugfindRDataset", core = core, sourceClass = sourceClass)
}
