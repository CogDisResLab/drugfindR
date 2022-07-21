# library(readr)
# library(dplyr)
# library(tibble)
# library(purrr)
# library(tidyr)
#
# empty_signature <- tibble::tibble(
#   signatureID = rep(NA, 978),
#   ID_geneid = rep(NA, 978),
#   Name_GeneSymbol = rep(NA, 978),
#   Value_LogDiffExp = rep(NA, 978),
#   Significance_pvalue = rep(NA, 978)
# )
#
# col_names <- colnames(empty_signature)
#
# col_spec <- cols(
#   signatureID = col_character(),
#   ID_geneid = col_double(),
#   Name_GeneSymbol = col_character(),
#   Value_LogDiffExp = col_double(),
#   Significance_pvalue = col_double()
# )
#
# kd_signature <- readr::read_csv("reference/LINCSKD_28.csv", col_types = col_spec)[]
# cp_signature <- readr::read_csv("reference/LINCSCP_5821.csv", col_types = col_spec)[]
# oe_signature <- readr::read_csv("reference/LINCSOE_104.csv", col_types = col_spec)[]
#
# test_that("empty signature when filtered is empty", {
#   expect_equal(nrow(filter_signature(empty_signature)), 0)
# })
#
# # Checking filtering of knockdown signatures
#
# test_that("knockdown signature has correct rows when filtering in any direction", {
#   expect_equal(nrow(filter_signature(kd_signature, "any", 0)), nrow(read_csv("reference/kd_signature_any_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "any", 0.5)), nrow(read_csv("reference/kd_signature_any_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "any", 0.85)), nrow(read_csv("reference/kd_signature_any_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "any", 1)), nrow(read_csv("reference/kd_signature_any_1.csv", col_types = col_spec)[]))
# })
#
# test_that("knockdown signature has correct rows when filtering in up direction", {
#   expect_equal(nrow(filter_signature(kd_signature, "up", 0)), nrow(read_csv("reference/kd_signature_up_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "up", 0.5)), nrow(read_csv("reference/kd_signature_up_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "up", 0.85)), nrow(read_csv("reference/kd_signature_up_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "up", 1)), nrow(read_csv("reference/kd_signature_up_1.csv", col_types = col_spec)[]))
# })
#
# test_that("knockdown signature has correct rows when filtering in down direction", {
#   expect_equal(nrow(filter_signature(kd_signature, "down", 0)), nrow(read_csv("reference/kd_signature_down_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "down", 0.5)), nrow(read_csv("reference/kd_signature_down_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "down", 0.85)), nrow(read_csv("reference/kd_signature_down_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(kd_signature, "down", 1)), nrow(read_csv("reference/kd_signature_down_1.csv", col_types = col_spec)[]))
# })
#
# test_that("knockdown filtered signatures have the correct columns", {
#   expect_equal(names(filter_signature(kd_signature, "any", 0)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "any", 0.5)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "any", 0.85)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "any", 1)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "up", 0)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "up", 0.5)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "up", 0.85)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "up", 1)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "down", 0)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "down", 0.5)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "down", 0.85)), col_names)
#   expect_equal(names(filter_signature(kd_signature, "down", 1)), col_names)
# })
#
# test_that("knockdown signature has correct content when filtering in any direction", {
#   expect_equal(filter_signature(kd_signature, "any", 0), read_csv("reference/kd_signature_any_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "any", 0.5), read_csv("reference/kd_signature_any_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "any", 0.85), read_csv("reference/kd_signature_any_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "any", 1), read_csv("reference/kd_signature_any_1.csv", col_types = col_spec)[])
# })
#
# test_that("knockdown signature has correct content when filtering in up direction", {
#   expect_equal(filter_signature(kd_signature, "up", 0), read_csv("reference/kd_signature_up_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "up", 0.5), read_csv("reference/kd_signature_up_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "up", 0.85), read_csv("reference/kd_signature_up_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "up", 1), read_csv("reference/kd_signature_up_1.csv", col_types = col_spec)[])
# })
#
# test_that("knockdown signature has correct content when filtering in down direction", {
#   expect_equal(filter_signature(kd_signature, "down", 0), read_csv("reference/kd_signature_down_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "down", 0.5), read_csv("reference/kd_signature_down_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "down", 0.85), read_csv("reference/kd_signature_down_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(kd_signature, "down", 1), read_csv("reference/kd_signature_down_1.csv", col_types = col_spec)[])
# })
#
# # Checking filtering of overexpression signatures
#
# test_that("overexpression signature has correct rows when filtering in any direction", {
#   expect_equal(nrow(filter_signature(oe_signature, "any", 0)), nrow(read_csv("reference/oe_signature_any_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "any", 0.5)), nrow(read_csv("reference/oe_signature_any_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "any", 0.85)), nrow(read_csv("reference/oe_signature_any_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "any", 1)), nrow(read_csv("reference/oe_signature_any_1.csv", col_types = col_spec)[]))
# })
#
# test_that("overexpression signature has correct rows when filtering in up direction", {
#   expect_equal(nrow(filter_signature(oe_signature, "up", 0)), nrow(read_csv("reference/oe_signature_up_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "up", 0.5)), nrow(read_csv("reference/oe_signature_up_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "up", 0.85)), nrow(read_csv("reference/oe_signature_up_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "up", 1)), nrow(read_csv("reference/oe_signature_up_1.csv", col_types = col_spec)[]))
# })
#
# test_that("overexpression signature has correct rows when filtering in down direction", {
#   expect_equal(nrow(filter_signature(oe_signature, "down", 0)), nrow(read_csv("reference/oe_signature_down_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "down", 0.5)), nrow(read_csv("reference/oe_signature_down_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "down", 0.85)), nrow(read_csv("reference/oe_signature_down_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(oe_signature, "down", 1)), nrow(read_csv("reference/oe_signature_down_1.csv", col_types = col_spec)[]))
# })
#
# test_that("overexpression filtered signatures have the correct columns", {
#   expect_equal(names(filter_signature(oe_signature, "any", 0)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "any", 0.5)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "any", 0.85)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "any", 1)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "up", 0)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "up", 0.5)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "up", 0.85)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "up", 1)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "down", 0)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "down", 0.5)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "down", 0.85)), col_names)
#   expect_equal(names(filter_signature(oe_signature, "down", 1)), col_names)
# })
#
# test_that("overexpression signature has correct content when filtering in any direction", {
#   expect_equal(filter_signature(oe_signature, "any", 0), read_csv("reference/oe_signature_any_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "any", 0.5), read_csv("reference/oe_signature_any_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "any", 0.85), read_csv("reference/oe_signature_any_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "any", 1), read_csv("reference/oe_signature_any_1.csv", col_types = col_spec)[])
# })
#
# test_that("overexpression signature has correct content when filtering in up direction", {
#   expect_equal(filter_signature(oe_signature, "up", 0), read_csv("reference/oe_signature_up_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "up", 0.5), read_csv("reference/oe_signature_up_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "up", 0.85), read_csv("reference/oe_signature_up_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "up", 1), read_csv("reference/oe_signature_up_1.csv", col_types = col_spec)[])
# })
#
# test_that("overexpression signature has correct content when filtering in down direction", {
#   expect_equal(filter_signature(oe_signature, "down", 0), read_csv("reference/oe_signature_down_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "down", 0.5), read_csv("reference/oe_signature_down_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "down", 0.85), read_csv("reference/oe_signature_down_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(oe_signature, "down", 1), read_csv("reference/oe_signature_down_1.csv", col_types = col_spec)[])
# })
#
# # Checking filtering of chemical perturbagen signatures
#
# test_that("chemical perturbagen signature has correct rows when filtering in any direction", {
#   expect_equal(nrow(filter_signature(cp_signature, "any", 0)), nrow(read_csv("reference/cp_signature_any_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "any", 0.5)), nrow(read_csv("reference/cp_signature_any_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "any", 0.85)), nrow(read_csv("reference/cp_signature_any_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "any", 1)), nrow(read_csv("reference/cp_signature_any_1.csv", col_types = col_spec)[]))
# })
#
# test_that("chemical perturbagen signature has correct rows when filtering in up direction", {
#   expect_equal(nrow(filter_signature(cp_signature, "up", 0)), nrow(read_csv("reference/cp_signature_up_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "up", 0.5)), nrow(read_csv("reference/cp_signature_up_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "up", 0.85)), nrow(read_csv("reference/cp_signature_up_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "up", 1)), nrow(read_csv("reference/cp_signature_up_1.csv", col_types = col_spec)[]))
# })
#
# test_that("chemical perturbagen signature has correct rows when filtering in down direction", {
#   expect_equal(nrow(filter_signature(cp_signature, "down", 0)), nrow(read_csv("reference/cp_signature_down_0.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "down", 0.5)), nrow(read_csv("reference/cp_signature_down_0.5.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "down", 0.85)), nrow(read_csv("reference/cp_signature_down_0.85.csv", col_types = col_spec)[]))
#   expect_equal(nrow(filter_signature(cp_signature, "down", 1)), nrow(read_csv("reference/cp_signature_down_1.csv", col_types = col_spec)[]))
# })
#
# test_that("chemical perturbagen filtered signatures have the correct columns", {
#   expect_equal(names(filter_signature(cp_signature, "any", 0)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "any", 0.5)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "any", 0.85)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "any", 1)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "up", 0)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "up", 0.5)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "up", 0.85)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "up", 1)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "down", 0)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "down", 0.5)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "down", 0.85)), col_names)
#   expect_equal(names(filter_signature(cp_signature, "down", 1)), col_names)
# })
#
# test_that("chemical perturbagen signature has correct content when filtering in any direction", {
#   expect_equal(filter_signature(cp_signature, "any", 0), read_csv("reference/cp_signature_any_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "any", 0.5), read_csv("reference/cp_signature_any_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "any", 0.85), read_csv("reference/cp_signature_any_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "any", 1), read_csv("reference/cp_signature_any_1.csv", col_types = col_spec)[])
# })
#
# test_that("chemical perturbagen signature has correct content when filtering in up direction", {
#   expect_equal(filter_signature(cp_signature, "up", 0), read_csv("reference/cp_signature_up_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "up", 0.5), read_csv("reference/cp_signature_up_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "up", 0.85), read_csv("reference/cp_signature_up_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "up", 1), read_csv("reference/cp_signature_up_1.csv", col_types = col_spec)[])
# })
#
# test_that("chemical perturbagen signature has correct content when filtering in down direction", {
#   expect_equal(filter_signature(cp_signature, "down", 0), read_csv("reference/cp_signature_down_0.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "down", 0.5), read_csv("reference/cp_signature_down_0.5.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "down", 0.85), read_csv("reference/cp_signature_down_0.85.csv", col_types = col_spec)[])
#   expect_equal(filter_signature(cp_signature, "down", 1), read_csv("reference/cp_signature_down_1.csv", col_types = col_spec)[])
# })
