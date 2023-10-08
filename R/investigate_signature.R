#' Investigate a given DGE dataset
#'
#' `r lifecycle::badge("experimental")`
#'
#' This function takes a DGE Data frame and then
#' finds concordant signatures to that.
#' This generates an L1000 signature from the DGE
#' dataset and then uploads that signature to
#' iLINCS to find the relevant concordant (or discordant) signatures
#'
#' @param expr A dataframe that has differential gene expression analysis
#' @param output_lib The library to search
#' @param filter_threshold The Filtering threshold.
#' @param filter_prop The Filtering proportion.
#' @param similarity_threshold The Similarity Threshold
#' @param paired Logical. Whether to query iLINCS separately
#' for up and down regulated genes
#' @param output_cell_lines A character vector of cell lines
#' to restrict the output search to.
#' @param gene_column The name of the column that has gene symbols
#' @param logfc_column The name of the column that has log_2 fold-change values
#' @param pval_column  The name of the column that has p-values
#' @param source_name (Optional) An annotation column to identify
#' the signature by name
#' @param source_cell_line (Optional) An annotation column to specify
#' the cell line for the input data
#' @param source_time (Optional) An annotation column to specify the
#' time for the input data
#' @param source_concentration (Optional) An annotation column to specify
#' the concentration for the input data
#'
#' @return A tibble with the the similarity scores and signature metadata
#' @export
#'
#' @importFrom dplyr mutate select any_of
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
investigate_signature <- function(expr,
                                  output_lib,
                                  filter_threshold = NULL,
                                  filter_prop = NULL,
                                  similarity_threshold = 0.2,
                                  paired = TRUE,
                                  output_cell_lines = NULL,
                                  gene_column = "Symbol",
                                  logfc_column = "logFC",
                                  pval_column = "PValue",
                                  source_name = "Input",
                                  source_cell_line = "NA",
                                  source_time = "NA",
                                  source_concentration = "NA") {
    libs <- c("OE", "KD", "CP")

    if (!output_lib %in% libs) {
        stop("Output library must be one of 'OE', 'KD', 'CP'")
    }

    if (missing(output_lib)) {
        stop("Please specify an output library")
    }

    expr_signature <- expr %>%
        prepare_signature(
            gene_column = gene_column,
            logfc_column = logfc_column,
            pval_column = pval_column
        )

    signature_id <- unique(expr_signature[["signatureID"]])

    if (paired) {
        filtered_up <- expr_signature %>%
            filter_signature(
                direction = "up",
                threshold = filter_threshold,
                prop = filter_prop
            )

        filtered_down <- expr_signature %>%
            filter_signature(
                direction = "down",
                threshold = filter_threshold,
                prop = filter_prop
            )

        concordant_up <- filtered_up %>%
            get_concordants(ilincs_library = output_lib, sig_direction = "Up")

        concordant_down <- filtered_down %>%
            get_concordants(ilincs_library = output_lib, sig_direction = "Down")


        consensus_targets <-
            consensus_concordants(
                concordant_up,
                concordant_down,
                paired = paired,
                cell_line = output_cell_lines,
                cutoff = similarity_threshold
            )
    } else {
        filtered <- expr_signature %>%
            filter_signature(
                direction = "any",
                threshold = filter_threshold,
                prop = filter_prop
            )

        concordants <- filtered %>%
            get_concordants(ilincs_library = output_lib)

        consensus_targets <-
            consensus_concordants(
                concordants,
                paired = paired,
                cell_line = output_cell_lines,
                cutoff = similarity_threshold
            )
    }

    augmented <- consensus_targets %>%
        dplyr::mutate(
            SourceSignature = signature_id,
            Source = source_name,
            SourceCellLine = source_cell_line,
            SourceTime = source_time,
            InputSignatureDirection = sig_direction
        ) %>%
        dplyr::select(
            dplyr::any_of(
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
            )
        )

    augmented
}
