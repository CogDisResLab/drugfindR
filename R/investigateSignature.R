#' @include utilities.R
#' @include getSignature.R prepareSignature.R
#' @include getConcordants.R consensusConcordants.R filterSignature.R
NULL

#' Investigate a given DGE dataset
#' `r lifecycle::badge("experimental")`
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
#'
#' # Investigate a signature
#'
#' # Load and prepare the signature
#' inputSignature <- read.table(system.file("extdata",
#'     "dCovid_diffexp.tsv",
#'     package = "drugfindR"
#' ), header = TRUE)
#'
#'
#' # Investigate the signature
#'
#' investigatedSignature <- investigateSignature(inputSignature,
#'     outputLib = "CP",
#'     filterThreshold = 0.5,
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC",
#'     pvalColumn = "PValue"
#' )
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
    sourceConcentration = "NA") {
    stopIfInvalidLibraries(outputLib)
    if (missing(outputLib)) {
        stop("Please specify an output library")
    }

    exprSignature <- expr %>%
        prepareSignature(
            geneColumn = geneColumn,
            logfcColumn = logfcColumn,
            pvalColumn = pvalColumn
        )

    signatureId <- unique(exprSignature[["signatureID"]])

    if (paired) {
        filteredUp <- exprSignature %>%
            filterSignature(
                direction = "up",
                threshold = filterThreshold,
                prop = filterProp
            )

        filteredDown <- exprSignature %>%
            filterSignature(
                direction = "down",
                threshold = filterThreshold,
                prop = filterProp
            )

        concordantUp <- filteredUp %>%
            getConcordants(ilincsLibrary = outputLib)

        concordantDown <- filteredDown %>%
            getConcordants(ilincsLibrary = outputLib)


        consensusTargets <-
            consensusConcordants(
                concordantUp,
                concordantDown,
                paired = paired,
                cellLine = outputCellLines,
                cutoff = similarityThreshold
            )
    } else {
        filtered <- exprSignature %>%
            filterSignature(
                direction = "any",
                threshold = filterThreshold,
                prop = filterProp
            )

        concordants <- filtered %>%
            getConcordants(ilincsLibrary = outputLib)

        consensusTargets <-
            consensusConcordants(
                concordants,
                paired = paired,
                cellLine = outputCellLines,
                cutoff = similarityThreshold
            )
    }

    augmented <- consensusTargets %>%
        dplyr::mutate(
            SourceSignature = signatureId,
            Source = sourceName,
            SourceCellLine = sourceCellLine,
            SourceTime = sourceTime,
        ) %>%
        dplyr::select(
            dplyr::any_of(c(
                "Source",
                "Target",
                "Similarity",
                "SourceSignature",
                "SourceCellLine",
                "InputSignatureDirection",
                "SourceConcentration",
                "SourceTime",
                "TargetSignature",
                "TargetCellLine",
                "TargetConcentration",
                "TargetTime"
            ))
        )

    augmented
}
