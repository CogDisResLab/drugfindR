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
    rds_path <- file.path(test_path(), "fixtures", "example_signature.RDS")
    if (file.exists(rds_path)) {
        readr::read_rds(rds_path)
    } else {
        get_signature("LINCSKD_28") |>
            saveRDS(file = rds_path)
    }
}

## Generate concordants for a signature

concordants_cp <- function() {
    rds_path <- file.path(test_path(), "fixtures", "concordants_cp.RDS")
    if (file.exists(rds_path)) {
        readr::read_rds(rds_path)
    } else {
        get_concordants(
            {
                example_signature() |> filter_signature(threshold = 1.0)
            },
            "CP",
            "any"
        ) |>
            saveRDS(file = rds_path, compress = "xz")
    }
}

concordants_cp_paired <- function() {
    rds_path <- file.path(test_path(), "fixtures", "concordants_cp_paired.RDS")
    if (file.exists(rds_path)) {
        readr::read_rds(rds_path)
    } else {
        signature_upregulated <- example_signature() |> filter_signature(threshold = 1.0, direction = "up")
        signature_downregulated <- example_signature() |> filter_signature(threshold = 1.0, direction = "down")
        up_concordants <- get_concordants(
            signature_upregulated,
            "CP",
            "up"
        )
        down_concordants <- get_concordants(
            signature_downregulated,
            "CP",
            "down"
        )
        list(up_concordants, down_concordants) |>
            saveRDS(file = rds_path, compress = "xz")
    }
}

concordants_oe <- function() {
    rds_path <- file.path(test_path(), "fixtures", "concordants_oe.RDS")
    if (file.exists(rds_path)) {
        readr::read_rds(rds_path)
    } else {
        get_concordants(
            {
                example_signature() |> filter_signature(threshold = 1.0)
            },
            "OE",
            "any"
        ) |>
            saveRDS(file = rds_path, compress = "xz")
    }
}

# Concordants Column Names

concordants_col_names <- function() {
    colnames(concordants_cp())
}

## Consensus CP Concordants Column Names

consensus_concordants_col_names <- function() { # nolint: object_length_linter.
    colnames(consensus_concordants(concordants_cp()))
}

## Consensus OE Concordants Column Names

consensus_concordants_oe_col_names <- function() { # nolint: object_length_linter.
    colnames(consensus_concordants(concordants_oe()))
}
