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
#'   1. Ensures all column names are character strings
#'   1. Validates that specified columns exist in the input dataframe
#'   1. Checks that the dataframe is not empty
#'
#' @keywords internal
#'
#' @examples NULL
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
#' @importFrom dplyr inner_join mutate select any_of
#' @importFrom rlang .data
#' @import DFplyr
.mapToL1000WithPvalues <- function(filteredData, geneColumn, logfcColumn, pvalColumn) {
    # Convert DataFrame to tibble for dplyr operations
    # This avoids scoping issues with all_of() in DFplyr
    filteredData <- tibble::as_tibble(filteredData)

    # Select only the relevant columns from the input to avoid
    # carrying through any pre-existing L1000-style columns that
    # would collide with our target names after renaming.
    filteredSubset <- filteredData |>
        dplyr::select(dplyr::all_of(c(geneColumn, logfcColumn, pvalColumn)))

    l1000 |> # nolint: object_usage_linter.
        dplyr::inner_join(filteredSubset,
            by = c(SYMBOL = geneColumn), relationship = "many-to-many"
        ) |>
        # NOTE: For the renaming, we are using an in-place
        # variable by unquoting a literal string.
        # This avoids the "No global item" note
        # in `R CMD CHECK`.
        dplyr::rename(
            ID_geneid = !!"ENTREZID",
            Name_GeneSymbol = !!"L1000",
            Value_LogDiffExp = !!logfcColumn,
            Significance_pvalue = !!pvalColumn
        ) |>
        dplyr::mutate(signatureID = "InputSig") |>
        dplyr::select(dplyr::any_of(c(
            "signatureID",
            "ID_geneid",
            "Name_GeneSymbol",
            "Value_LogDiffExp",
            "Significance_pvalue"
        ))) |>
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
#' @importFrom dplyr inner_join mutate select any_of
#' @importFrom rlang .data
#' @import DFplyr
.mapToL1000WithoutPvalues <- function(filteredData, geneColumn, logfcColumn) {
    # Convert DataFrame to tibble for dplyr operations
    # This avoids scoping issues with all_of() in DFplyr
    filteredData <- tibble::as_tibble(filteredData)

    # Select only the relevant columns from the input to avoid
    # carrying through any pre-existing L1000-style columns that
    # would collide with our target names after renaming.
    filteredSubset <- filteredData |>
        dplyr::select(dplyr::all_of(c(geneColumn, logfcColumn)))

    l1000 |> # nolint: object_usage_linter.
        dplyr::inner_join(filteredSubset,
            by = c(SYMBOL = geneColumn), relationship = "many-to-many"
        ) |>
        # NOTE: For the renaming, we are using an in-place
        # variable by unquoting a literal string.
        # This avoids the "No global item" note
        # in `R CMD CHECK`.
        dplyr::rename(
            ID_geneid = !!"ENTREZID",
            Name_GeneSymbol = !!"L1000",
            Value_LogDiffExp = !!logfcColumn
        ) |>
        dplyr::mutate(signatureID = "InputSig") |>
        dplyr::select(dplyr::any_of(c(
            "signatureID",
            "ID_geneid",
            "Name_GeneSymbol",
            "Value_LogDiffExp",
            "Significance_pvalue"
        ))) |>
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
#'   1. `.mapToL1000WithPvalues` when p-value column is specified
#'   1. `.mapToL1000WithoutPvalues` when p-value column is NA
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
#' @importFrom dplyr filter select any_of inner_join mutate
#' @importFrom rlang .data
#'
#' @examples
#' # Load example differential expression data from package
#' dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
#'     package = "drugfindR"
#' )
#' dge_data <- read.delim(dge_file)
#'
#' # Prepare signature with p-values (standard workflow)
#' signature <- prepareSignature(
#'     dge_data,
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC",
#'     pvalColumn = "PValue"
#' )
#' head(signature)
#'
#' # Prepare signature without p-values
#' signature_no_pval <- prepareSignature(
#'     dge_data,
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC",
#'     pvalColumn = NA
#' )
#' head(signature_no_pval)
#'
#' # Custom column names example
#' custom_dge <- data.frame(
#'     Gene = c("TP53", "MYC", "BRCA1", "EGFR"),
#'     FC = c(2.5, -1.8, 3.2, -2.1),
#'     Pval = c(0.001, 0.01, 0.0001, 0.005)
#' )
#'
#' custom_signature <- prepareSignature(
#'     custom_dge,
#'     geneColumn = "Gene",
#'     logfcColumn = "FC",
#'     pvalColumn = "Pval"
#' )
#' print(custom_signature)
prepareSignature <- function(
    dge, geneColumn = "Symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue") {
    # Validate input parameters
    .validatePrepareSignatureInput(dge, geneColumn, logfcColumn, pvalColumn)

    # Filter data to L1000 genes and process into signature format
    .processToL1000Signature(dge, geneColumn, logfcColumn, pvalColumn)
}
