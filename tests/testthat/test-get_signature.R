library(readr)
library(dplyr)
library(tibble)
library(purrr)

empty_signature <- tibble::tibble(
  signatureID = rep(NA, 978),
  ID_geneid = rep(NA, 978),
  Name_GeneSymbol = rep(NA, 978),
  Value_LogDiffExp = rep(NA, 978),
  Significance_pvalue = rep(NA, 978)
)

col_names <- colnames(empty_signature)

col_spec <- cols(
  signatureID = col_character(),
  ID_geneid = col_double(),
  Name_GeneSymbol = col_character(),
  Value_LogDiffExp = col_double(),
  Significance_pvalue = col_double()
)

kd_signature_known <- readr::read_csv("reference/LINCSKD_28.csv", col_types = col_spec)[]
cp_signature_known <- readr::read_csv("reference/LINCSCP_5821.csv", col_types = col_spec)[]
oe_signature_known <- readr::read_csv("reference/LINCSOE_104.csv", col_types = col_spec)[]

inv_signature <- get_signature("LINCS_INV")
kd_signature <- get_signature("LINCSKD_28")
cp_signature <- get_signature("LINCSCP_5821")
oe_signature <- get_signature("LINCSOE_104")

# Testing invalid signature
test_that("correct number of rows for the invalid signature", {
  expect_equal(nrow(inv_signature), 978)
})

test_that("everything NA for invalid signature", {
  expect_equal(
    all(purrr::flatten_lgl(purrr::map(inv_signature, is.na))),
    TRUE)
})

test_that("correct columns for the invalid signature", {
  expect_equal(names(inv_signature), col_names)
})

test_that("correct content for the invalid signature", {
  expect_equal(inv_signature, empty_signature)
})

# Testing the knockdown signature
test_that("correct number of rows for the knockdown signature", {
  expect_equal(nrow(kd_signature), 978)
})

test_that("correct columns for the knockdown signature", {
  expect_equal(names(kd_signature), col_names)
})

test_that("nothing NA for knockdown signature", {
  expect_equal(
    any(purrr::flatten_lgl(purrr::map(kd_signature, is.na))),
    FALSE)
})

test_that("correct content for the knockdown signature", {
  expect_equal(kd_signature, kd_signature_known)
})

# Testing the chemical perturbagen signature
test_that("correct number of rows for the chemical perturbagen signature", {
  expect_equal(nrow(cp_signature), 978)
})

test_that("correct columns for the chemical perturbagen signature", {
  expect_equal(names(cp_signature), col_names)
})

test_that("nothing NA for chemical perturbagen signature", {
  expect_equal(
    any(purrr::flatten_lgl(purrr::map(cp_signature, is.na))),
    FALSE)
})

test_that("correct content for the chemical perturbagen signature", {
  expect_equal(cp_signature, cp_signature_known)
})

# Testing the Overexpression signature
test_that("correct number of rows for the overexpression signature", {
  expect_equal(nrow(oe_signature), 978)
})

test_that("correct columns for the overexpression signature", {
  expect_equal(names(oe_signature), col_names)
})

test_that("nothing NA for overexpression signature", {
  expect_equal(
    any(purrr::flatten_lgl(purrr::map(oe_signature, is.na))),
    FALSE)
})

test_that("correct content for the overexpression signature", {
  expect_equal(oe_signature, oe_signature_known)
})
