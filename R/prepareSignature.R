#' Validate prepareSignature input parameters
#'
#' This internal function validates all input parameters for the prepareSignature
#' function to ensure they meet the required constraints.
#'
#' @param dge A dataframe-like object containing differential gene expression data.
#' @param geneColumn Character string specifying the column name containing gene symbols.
#' @param logfcColumn Character string specifying the column name containing log fold-change values.
#' @param pvalColumn Character string specifying the column name containing p-values, or NA.
#'
#' @return Invisible NULL. The function throws an error if validation fails.
#'
#' @details
#' This function performs the following validations:
#' \enumerate{
#'   \item Ensures all column names are character strings
#'   \item Validates that specified columns exist in the input dataframe
#'   \item Checks that the dataframe is not empty
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Valid calls (no errors)
#' data <- data.frame(Symbol = c("TP53", "MYC"), logFC = c(1.5, -2.1), PValue = c(0.01, 0.05))
#' .validatePrepareSignatureInput(data, "Symbol", "logFC", "PValue")
#' .validatePrepareSignatureInput(data, "Symbol", "logFC", NA)
#'
#' # Invalid calls (will throw errors)
#' .validatePrepareSignatureInput(data, "NonExistent", "logFC", "PValue")
#' .validatePrepareSignatureInput(data, "Symbol", "NonExistent", "PValue")
#' }
.validatePrepareSignatureInput <- function(dge, geneColumn, logfcColumn, pvalColumn) {
    # Validate input types
    if (!is.character(geneColumn) || length(geneColumn) != 1L) {
        stop("geneColumn must be a single character string", call. = FALSE)
    }

    if (!is.character(logfcColumn) || length(logfcColumn) != 1L) {
        stop("logfcColumn must be a single character string", call. = FALSE)
    }

    if (!is.character(pvalColumn) && !is.na(pvalColumn)) {
        stop("pvalColumn must be a character string or NA", call. = FALSE)
    }

    if (is.character(pvalColumn) && length(pvalColumn) != 1L) {
        stop("pvalColumn must be a single character string or NA", call. = FALSE)
    }

    # Validate dataframe
    if (is.null(dge) || nrow(dge) == 0L) {
        stop("dge must be a non-empty dataframe-like object", call. = FALSE)
    }

    # Validate column existence
    if (!geneColumn %in% names(dge)) {
        stop("geneColumn '", geneColumn, "' not found in the dataframe", call. = FALSE)
    }

    if (!logfcColumn %in% names(dge)) {
        stop("logfcColumn '", logfcColumn, "' not found in the dataframe", call. = FALSE)
    }

    if (!is.na(pvalColumn) && !pvalColumn %in% names(dge)) {
        stop("pvalColumn '", pvalColumn, "' not found in the dataframe", call. = FALSE)
    }
}

#' Filter differential expression data to L1000 genes
#'
#' This internal function filters the input differential expression data to include
#' only genes present in the L1000 gene set.
#'
#' @param dge A dataframe containing differential gene expression data.
#' @param geneColumn Character string specifying the column name containing gene symbols.
#' @param logfcColumn Character string specifying the column name containing log fold-change values.
#' @param pvalColumn Character string specifying the column name containing p-values, or NA.
#'
#' @return A filtered dataframe containing only L1000 genes with the specified columns.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Filters the input data to genes present in the L1000 gene set
#'   \item Selects only the required columns (gene, logFC, and optionally p-value)
#'   \item Returns the filtered dataset for further processing
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr filter select any_of
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' data <- data.frame(Symbol = c("TP53", "MYC", "INVALID"), logFC = c(1.5, -2.1, 0.5))
#' filtered <- .filterToL1000Genes(data, "Symbol", "logFC", NA)
#' }
.filterToL1000Genes <- function(dge, geneColumn, logfcColumn, pvalColumn = NA) {
    if (!is.na(pvalColumn)) {
        dge %>%
            dplyr::filter(.data[[geneColumn]] %in% l1000[["SYMBOL"]]) %>%
            dplyr::select(dplyr::any_of(c(geneColumn, logfcColumn, pvalColumn)))
    } else {
        dge %>%
            dplyr::filter(.data[[geneColumn]] %in% l1000[["SYMBOL"]]) %>%
            dplyr::select(dplyr::any_of(c(geneColumn, logfcColumn)))
    }
}

#' Map filtered data to L1000 format with p-values
#'
#' This internal function maps the filtered differential expression data to the
#' standardized L1000 signature format, including p-value information.
#'
#' @param filteredData A dataframe containing filtered differential expression data.
#' @param geneColumn Character string specifying the column name containing gene symbols.
#' @param logfcColumn Character string specifying the column name containing log fold-change values.
#' @param pvalColumn Character string specifying the column name containing p-values.
#'
#' @return A tibble with standardized L1000 signature format including p-values.
#'
#' @keywords internal
#'
#' @importFrom dplyr inner_join rename mutate select any_of
#' @importFrom rlang .data
.mapToL1000WithPvalues <- function(filteredData, geneColumn, logfcColumn, pvalColumn) {
    l1000 %>%
        dplyr::inner_join(filteredData, by = c(SYMBOL = geneColumn)) %>%
        dplyr::rename(
            ID_geneid = !!"ENTREZID",
            Name_GeneSymbol = !!"L1000",
            Value_LogDiffExp = !!logfcColumn,
            Significance_pvalue = !!pvalColumn
        ) %>%
        dplyr::mutate(signatureID = "InputSig") %>%
        dplyr::select(dplyr::any_of(c(
            "signatureID",
            "ID_geneid",
            "Name_GeneSymbol",
            "Value_LogDiffExp",
            "Significance_pvalue"
        ))) %>%
        unique()
}

#' Map filtered data to L1000 format without p-values
#'
#' This internal function maps the filtered expression data to the
#' standardized L1000 signature format, without p-value information.
#'
#' @param filteredData A dataframe containing filtered differential expression data.
#' @param geneColumn Character string specifying the column name containing gene symbols.
#' @param logfcColumn Character string specifying the column name containing log fold-change values.
#'
#' @return A tibble with standardized L1000 signature format without p-values.
#'
#' @keywords internal
#'
#' @importFrom dplyr inner_join rename mutate select any_of
#' @importFrom rlang .data
.mapToL1000WithoutPvalues <- function(filteredData, geneColumn, logfcColumn) {
    l1000 %>%
        dplyr::inner_join(filteredData, by = c(SYMBOL = geneColumn)) %>%
        dplyr::rename(
            ID_geneid = !!"ENTREZID",
            Name_GeneSymbol = !!"L1000",
            Value_LogDiffExp = !!logfcColumn
        ) %>%
        dplyr::mutate(signatureID = "InputSig") %>%
        dplyr::select(dplyr::any_of(c(
            "signatureID",
            "ID_geneid",
            "Name_GeneSymbol",
            "Value_LogDiffExp",
            "Significance_pvalue"
        ))) %>%
        unique()
}

#' Process differential expression data into L1000 signature format
#'
#' This internal function orchestrates the conversion of filtered differential
#' expression data into the standardized L1000 signature format.
#'
#' @param filteredData A dataframe containing filtered differential expression data.
#' @param geneColumn Character string specifying the column name containing gene symbols.
#' @param logfcColumn Character string specifying the column name containing log fold-change values.
#' @param pvalColumn Character string specifying the column name containing p-values, or NA.
#'
#' @return A tibble with the standardized L1000 signature format.
#'
#' @details
#' This function dispatches to appropriate mapping functions based on whether
#' p-value information is available:
#' \enumerate{
#'   \item \code{.mapToL1000WithPvalues} when p-value column is specified
#'   \item \code{.mapToL1000WithoutPvalues} when p-value column is NA
#' }
#'
#' @keywords internal
.processToL1000Signature <- function(filteredData, geneColumn, logfcColumn, pvalColumn = NA) {
    if (!is.na(pvalColumn)) {
        .mapToL1000WithPvalues(filteredData, geneColumn, logfcColumn, pvalColumn)
    } else {
        .mapToL1000WithoutPvalues(filteredData, geneColumn, logfcColumn)
    }
}

#' Prepare an L1000 Signature from a given differential gene expression output
#' `r lifecycle::badge("stable")`
#'
#' This function takes a differential gene expression output from any pipeline
#' like edgeR or DeSeq2 or any that give you the gene symbol,
#' log_2 fold-change and p-value
#' and transforms that into an L1000 signature for later processing.
#'
#' @param dge A dataframe-like object that has the differential
#' gene expression information
#' @param geneColumn The name of the column that has gene symbols
#' @param logfcColumn The name of the column that has log_2 fold-change values
#' @param pvalColumn  The name of the column that has p-values
#'
#' @return A tibble with the L1000 signature.
#' @export
#'
#' @importFrom dplyr filter select any_of inner_join rename mutate
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
#' # Prepare an L1000 signature from a differential gene expression output
#'
#' inputSignature <- read.table(system.file("extdata",
#'     "dCovid_diffexp.tsv",
#'     package = "drugfindR"
#' ), header = TRUE)
#'
#' signature <- prepareSignature(inputSignature,
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC", pvalColumn = "PValue"
#' )
#'
#' head(signature)
prepareSignature <- function(
    dge,
    geneColumn = "Symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue") {
    # Validate input parameters
    .validatePrepareSignatureInput(dge, geneColumn, logfcColumn, pvalColumn)

    # Filter data to L1000 genes and process into signature format
    .filterToL1000Genes(dge, geneColumn, logfcColumn, pvalColumn) %>%
        .processToL1000Signature(geneColumn, logfcColumn, pvalColumn)
}
