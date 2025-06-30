#' @include class-drugfindRResult.R
NULL

#' Construct a drugfindRResult Object
#'
#' This function wraps a `tibble`, `data.frame`, or `DataFrame` into a
#' `drugfindRResult` S4 object, for standardized manipulation in the
#' `drugfindR` pipeline.
#'
#' @param input A data-like object — `data.frame`, `tibble`, or `DataFrame`
#'
#' @return An S4 object of class `drugfindRResult`
#' @export
#'
#' @importFrom methods new
#'
#' @examples
#' df <- tibble::tibble(
#'     Gene = letters[1:5],
#'     Value_LogDiffExp = c(-1, -0.5, 0, 0.5, 1)
#' )
#' sig <- drugfindRResult(df)
#' sig
drugfindRResult <- function(dfdataset) {
    methods::new(
        "drugfindRResult",
        core = dfdataset@core,
        sourceClass = dfdataset@sourceClass
    )
}
