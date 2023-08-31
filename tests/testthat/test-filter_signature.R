library(readr)
library(dplyr)
library(tibble)
library(purrr)


# Setup

empty_signature <- tibble::tibble(
    signatureID = rep(NA, 978L),
    ID_geneid = rep(NA, 978L),
    Name_GeneSymbol = rep(NA, 978L),
    Value_LogDiffExp = rep(NA, 978L),
    Significance_pvalue = rep(NA, 978L)
)

col_names <- colnames(empty_signature)

col_spec <- cols(
    signatureID = col_character(),
    ID_geneid = col_double(),
    Name_GeneSymbol = col_character(),
    Value_LogDiffExp = col_double(),
    Significance_pvalue = col_double()
)

signature_id <- "LINCSKD_28"

sig_dir <- tempdir()
sig_file <-
    file.path(sig_dir, stringr::str_glue("{signature_id}.csv"))

signature <-
    readr::read_csv(
        stringr::str_glue("reference/{signature_id}.csv"),
        col_types = col_spec
    )



# Test Input Validation

test_that("Not specifying threshold and prop causes error", {
    expect_error(filter_signature(kd_signature))
})

test_that("Specifying both threshold and prop causes error", {
    expect_error(filter_signature(kd_signature, threshold = 0.1, prop = 0.1))
})

test_that("empty signature when filtered is empty", {
    expect_identical(nrow(filter_signature(empty_signature, threshold = 0.0)), 0L)
})

# Testing Filtering
