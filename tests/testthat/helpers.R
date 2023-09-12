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
        readr::read_csv(file.path(test_path(), "fixtures", "example_signature.csv"), col_types = "cccnn")
    } else {
        get_signature("LINCSKD_28") |>
            readr::write_csv(file.path(test_path(), "fixtures", "example_signature.csv"))
    }
}

## Generate concordants for a signature

concordants_cp <- function() {
    if (file.exists(file.path(test_path(), "fixtures", "concordants_cp.csv"))) {
        readr::read_csv(file.path(test_path(), "fixtures", "concordants_cp.csv"), col_types = "cccccnnc")
    } else {
        get_concordants(
            {
                example_signature() |> filter_signature(threshold = 1.0)
            },
            "CP",
            "any"
        ) |>
            readr::write_csv(file.path(test_path(), "fixtures", "concordants_cp.csv"))
    }
}
