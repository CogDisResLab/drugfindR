# library(dplyr)
# library(readr)
# library(tibble)
# library(tidyr)
#
# input_col_spec <- cols(
#   Gene_ID = col_character(),
#   Symbol = col_character(),
#   logFC = col_double(),
#   logCPM = col_double(),
#   F = col_double(),
#   PValue = col_double()
# )
#
# output_col_spec <- cols(
#   signatureID = col_character(),
#   ID_geneid = col_character(),
#   Name_GeneSymbol = col_character(),
#   Value_LogDiffExp = col_double(),
#   Significance_pvalue = col_double()
# )
#
# out_cols <- names(output_col_spec$cols)
#
# reference_input <- read_csv("reference/reference_rnaseq.csv", col_types = input_col_spec)
# reference_output <- read_csv("reference/reference_prepared_signature.csv", col_types = output_col_spec)
#
# test_that("Signature is in the proper format", {
#   example <- prepare_signature(reference_input, "Symbol", "logFC", "PValue")
#   expect_equal(nrow(example), 978)
#   expect_equal(example, reference_output, ignore_attr = TRUE)
#   })
#
# test_that("Correct Error messages are generated",  {
#   expect_error(prepare_signature(reference_input, gene_column = "Gene"), "gene_column should be present in the dataframe")
#   expect_error(prepare_signature(reference_input, logfc_column = "log2fc"), "logfc_column should be present in the dataframe")
#   expect_error(prepare_signature(reference_input, pval_column = "Gene"), "pval_column should be present in the dataframe")
# })
