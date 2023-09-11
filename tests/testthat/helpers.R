# Helper functions for testing

# Load packages
library(tibble)

## Create Empty Signature

empty_signature <- function() {
    tibble::tibble(
        signatureID = rep(NA, 978L),
        ID_geneid = rep(NA, 978L),
        Name_GeneSymbol = rep(NA, 978L),
        Value_LogDiffExp = rep(NA, 978L),
        Significance_pvalue = rep(NA, 978L)
    )
}

## Signature Column names

signature_col_names <- function() {
    colnames(empty_signature())
}

## Return an example signature

example_signature <- function() {
    if (file.exists(file.path(test_path(), "fixtures", "example_signature.csv"))) {
        readr::read_csv(file.path(test_path(), "fixtures", "example_signature.csv"))
    } else {
        get_signature("LINCSKD_28")
    }
}
