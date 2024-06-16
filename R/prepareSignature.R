#' Prepare an L1000 Signature froma given differential gene expression output
#' `r lifecycle::badge("experimental")`
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
#'     "dCovid_diffexp.tsv.tar.xz",
#'     package = "drugfindR"
#' ), header = TRUE)
#'
#' signature <- prepareSignature(inputSignature,
#'     geneColumn = "hgnc_symbol",
#'     logfcColumn = "logFC", pvalColumn = "PValue"
#' )
#'
#' head(signature)
prepareSignature <- function(dge,
                             geneColumn = "Symbol",
                             logfcColumn = "logFC",
                             pvalColumn = "PValue") {
    if (!geneColumn %in% names(dge)) {
        stop("geneColumn should be present in the dataframe")
    }

    if (!logfcColumn %in% names(dge)) {
        stop("logfcColumn should be present in the dataframe")
    }

    if (!pvalColumn %in% names(dge) && !is.na(pvalColumn)) {
        stop("pvalColumn should be present in the dataframe")
    }

    if (!is.na(pvalColumn)) {
        filteredL1000 <- dge %>%
            dplyr::filter(.data[[geneColumn]] %in% l1000[["SYMBOL"]]) %>%
            dplyr::select(
                dplyr::any_of(
                    c(geneColumn, logfcColumn, pvalColumn)
                )
            )

        signature <- l1000 %>%
            dplyr::inner_join(filteredL1000, by = c(SYMBOL = geneColumn)) %>%
            dplyr::rename(
                ID_geneid = !!"ENTREZID",
                Name_GeneSymbol = !!"L1000",
                Value_LogDiffExp = !!logfcColumn,
                Significance_pvalue = !!pvalColumn
            ) %>%
            dplyr::mutate(signatureID = "InputSig") %>%
            dplyr::select(
                any_of(c(
                    "signatureID",
                    "ID_geneid",
                    "Name_GeneSymbol",
                    "Value_LogDiffExp",
                    "Significance_pvalue"
                ))
            ) %>%
            unique()
    } else {
        filteredL1000 <- dge %>%
            dplyr::filter(!!geneColumn %in% l1000[["SYMBOL"]]) %>%
            dplyr::select(dplyr::any_of(c(geneColumn, logfcColumn)))

        signature <- l1000 %>%
            dplyr::inner_join(filteredL1000, by = c(SYMBOL = geneColumn)) %>%
            dplyr::rename(
                ID_geneid = !!"ENTREZID",
                Name_GeneSymbol = !!"L1000",
                Value_LogDiffExp = !!logfcColumn
            ) %>%
            dplyr::mutate(signatureID = "InputSig") %>%
            dplyr::select(
                any_of(c(
                    "signatureID",
                    "ID_geneid",
                    "Name_GeneSymbol",
                    "Value_LogDiffExp",
                    "Significance_pvalue"
                ))
            ) %>%
            unique()
    }


    signature
}
