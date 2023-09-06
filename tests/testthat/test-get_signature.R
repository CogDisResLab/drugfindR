# Get correct signature
#

library(readr)
library(stringr)
library(tibble)
library(purrr)

empty_signature <- tibble::tibble(
    signatureID = rep(NA, 978L),
    ID_geneid = rep(NA, 978L),
    Name_GeneSymbol = rep(NA, 978L),
    Value_LogDiffExp = rep(NA, 978L),
    Significance_pvalue = rep(NA, 978L)
)

col_names <- colnames(empty_signature)

col_spec <- cols(
    signatureID = readr::col_character(),
    ID_geneid = readr::col_double(),
    Name_GeneSymbol = readr::col_character(),
    Value_LogDiffExp = readr::col_double(),
    Significance_pvalue = readr::col_double()
)

kd_signature_id <- "LINCSKD_28"
cp_signature_id <- "LINCSCP_5821"
oe_signature_id <- "LINCSOE_104"

sig_dir <- tempdir()
kd_file <- file.path(sig_dir, stringr::str_glue("{kd_signature_id}.csv"))
cp_file <- file.path(sig_dir, stringr::str_glue("{cp_signature_id}.csv"))
oe_file <- file.path(sig_dir, stringr::str_glue("{oe_signature_id}.csv"))

inv_signature <- get_signature("LINCS_INV")
kd_signature <- get_signature(kd_signature_id) |>
    readr::write_csv(kd_file)
cp_signature <- get_signature(cp_signature_id) |>
    readr::write_csv(cp_file)
oe_signature <- get_signature(oe_signature_id) |>
    readr::write_csv(oe_file)


# Test the invalid signature
test_that("everything NA for invalid signature", {
    expect_true(
        all(purrr::flatten_lgl(purrr::map(inv_signature, is.na)))
    )
})

test_that("correct columns for the invalid signature", {
    expect_named(inv_signature, col_names)
})

test_that("correct content for the invalid signature", {
    expect_identical(inv_signature, empty_signature)
})

# Testing the knockdown signature
test_that("correct number of rows for the knockdown signature", {
    expect_identical(nrow(kd_signature), 978L)
})

test_that("correct columns for the knockdown signature", {
    expect_named(inv_signature, col_names)
})

test_that("nothing NA for knockdown signature", {
    expect_false(
        any(purrr::flatten_lgl(purrr::map(kd_signature, is.na)))
    )
})

test_that("correct content for the knockdown signature", {
    expect_snapshot_file(kd_file)
})

# Testing the chemical perturbagen signature
test_that("correct number of rows for the knockdown signature", {
    expect_identical(nrow(cp_signature), 978L)
})

test_that("correct columns for the knockdown signature", {
    expect_named(cp_signature, col_names)
})

test_that("nothing NA for knockdown signature", {
    expect_false(
        any(purrr::flatten_lgl(purrr::map(cp_signature, is.na)))
    )
})

test_that("correct content for the knockdown signature", {
    expect_snapshot_file(cp_file)
})

# Testing the overexpression signature
test_that("correct number of rows for the knockdown signature", {
    expect_identical(nrow(oe_signature), 978L)
})

test_that("correct columns for the knockdown signature", {
    expect_named(oe_signature, col_names)
})

test_that("nothing NA for knockdown signature", {
    expect_false(
        any(purrr::flatten_lgl(purrr::map(oe_signature, is.na)))
    )
})

test_that("correct content for the knockdown signature", {
    expect_snapshot_file(oe_file)
})
