#' @include utilities.R
#' @include getSignature.R prepareSignature.R
#' @include getConcordants.R consensusConcordants.R filterSignature.R
NULL

#' Investigate a given DGE dataset
#' `r lifecycle::badge("stable")`
#'
#' This function takes a DGE Data frame and then
#' finds concordant signatures to that.
#' This generates an L1000 signature from the DGE
#' dataset and then uploads that signature to
#' iLINCS to find the relevant concordant (or discordant) signatures
#'
#' @param expr A dataframe that has differential gene expression analysis
#' @param outputLib The library to search
#' @param filterThreshold The Filtering threshold.
#' @param filterProp The Filtering proportion.
#' @param similarityThreshold The Similarity Threshold
#' @param paired Logical. Whether to query iLINCS separately
#' for up and down regulated genes
#' @param outputCellLines A character vector of cell lines
#' to restrict the output search to.
#' @param geneColumn The name of the column that has gene symbols
#' @param logfcColumn The name of the column that has log_2 fold-change values
#' @param pvalColumn  The name of the column that has p-values
#' @param sourceName (Optional) An annotation column to identify
#' the signature by name
#' @param sourceCellLine (Optional) An annotation column to specify
#' the cell line for the input data
#' @param sourceTime (Optional) An annotation column to specify the
#' time for the input data
#' @param sourceConcentration (Optional) An annotation column to specify
#' the concentration for the input data
#'
#' @return A tibble with the the similarity scores and signature metadata
#' @export
#'
#' @importFrom dplyr mutate select any_of
#' @importFrom rlang .data
#'
#' @examples
#' # Input validation example (no API calls)
#' mockExpr <- data.frame(
#'     Symbol = c("TP53", "MYC"),
#'     logFC = c(2.5, -1.8),
#'     PValue = c(0.001, 0.01)
#' )
#'
#' # Validate library parameter (should produce error)
#' tryCatch(
#'     investigateSignature(mockExpr, outputLib = "INVALID"),
#'     error = function(e) message("Expected error: invalid library")
#' )
#'
#' \donttest{
#' # This function makes multiple API calls to iLINCS and may take several minutes
#'
#' # Load differential expression data
#' inputSignature <- read.table(
#'     system.file("extdata", "dCovid_diffexp.tsv", package = "drugfindR"),
#'     header = TRUE
#' )
#'
#' # Investigate the signature against chemical perturbagen library
#' investigatedSignature <- investigateSignature(
#'     inputSignature,
#'     outputLib = "CP",
#'     filterThreshold = 0.5,
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC",
#'     pvalColumn = "PValue"
#' )
#' head(investigatedSignature)
#' }
investigateSignature <- function(
  expr,
  outputLib,
  filterThreshold = NULL,
  filterProp = NULL,
  similarityThreshold = 0.2,
  paired = TRUE,
  outputCellLines = NULL,
  geneColumn = "Symbol",
  logfcColumn = "logFC",
  pvalColumn = "PValue",
  sourceName = "Input",
  sourceCellLine = "NA",
  sourceTime = "NA",
  sourceConcentration = "NA"
) {
    stopIfInvalidLibraries(outputLib) # nolint: object_usage_linter.
    if (missing(outputLib)) {
        stop("Please specify an output library", call. = FALSE)
    }

    # Prepare signature once; downstream functions handle validation & errors
    exprSignature <- prepareSignature( # nolint: object_usage_linter.
        expr,
        geneColumn = geneColumn,
        logfcColumn = logfcColumn,
        pvalColumn = pvalColumn
    )
    signatureId <- unique(exprSignature[["signatureID"]])

    consensus <- .computeConsensusFromSignature( # nolint: object_usage_linter.
        exprSignature,
        outputLib = outputLib,
        filterThreshold = filterThreshold,
        filterProp = filterProp,
        similarityThreshold = similarityThreshold,
        paired = paired,
        outputCellLines = outputCellLines
    )

    consensus |>
        dplyr::mutate(
            SourceSignature = signatureId,
            Source = sourceName,
            SourceCellLine = sourceCellLine,
            SourceTime = sourceTime,
            SourceConcentration = sourceConcentration
        ) |>
        dplyr::select(dplyr::any_of(c(
            "Source", "Target", "Similarity", "SourceSignature",
            "SourceCellLine", "SourceConcentration", "SourceTime",
            "TargetSignature", "TargetCellLine", "TargetConcentration", "TargetTime",
            "InputSigDirection", "SignatureType", "pValue"
        )))
}
