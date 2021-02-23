library(readr)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(stringr)

signature_col_spec <- cols(
  signatureID = col_character(),
  ID_geneid = col_double(),
  Name_GeneSymbol = col_character(),
  Value_LogDiffExp = col_double(),
  Significance_pvalue = col_double()
)

concordant_col_spec <- cols(
  .default = col_character(),
  similarity = col_double(),
  pValue = col_double()
)

input_filenames <- expand_grid(signature = c("kd_signature", "oe_signature", "cp_signature"), direction = c("any", "up", "down"), threshold = c(0.85)) %>% unite("filename") %>% pull(filename) %>% str_c(., "csv", sep = ".") %>% file.path("reference", .)

input_data <- map(input_filenames, ~ read_csv(.x, col_types = signature_col_spec))

reference_filenames <- expand_grid(concordant = "concordant", signature = c("kd_signature", "oe_signature", "cp_signature"), direction = c("any", "up", "down"), threshold = c(0.85), output_lib = c("OE", "KD", "CP")) %>% unite("filename") %>% pull(filename) %>% str_c(., "csv", sep = ".") %>% file.path("reference", .)

reference_data <- map(reference_filenames, ~ read_csv(.x, col_types = concordant_col_spec, na = c("")))

test_that("Error thrown at invalid signature", {
  expect_error(get_concordants("a"), "signature must be a data frame or data frame like object")
})

test_that("Error thrown at invalid library", {
  expect_error(get_concordants(tibble(), "AB"), "library must be one of 'OE', 'KD' or 'CP'")
  expect_error(get_concordants(tibble(), "CD"), "library must be one of 'OE', 'KD' or 'CP'")
  expect_error(get_concordants(tibble(), "OJ"), "library must be one of 'OE', 'KD' or 'CP'")
})

# Testing Knockdown Signatures
#
## Testing the results for the any filtered signature
test_that("Results from the up knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[1]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[1]]))
  expect_equal(r, reference_data[[1]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[1]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[2]]))
  expect_equal(r, reference_data[[2]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[1]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[3]]))
  expect_equal(r, reference_data[[3]], ignore_attr = TRUE)
})

## Testing the results for the up filtered signature
test_that("Results from the up knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[2]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[4]]))
  expect_equal(r, reference_data[[4]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[2]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[5]]))
  expect_equal(r, reference_data[[5]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[2]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[6]]))
  expect_equal(r, reference_data[[6]], ignore_attr = TRUE)
})

## Testing the results for the down filtered signature
test_that("Results from the down knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[3]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[7]]))
  expect_equal(r, reference_data[[7]], ignore_attr = TRUE)
})

test_that("Results from the down knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[3]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[8]]))
  expect_equal(r, reference_data[[8]], ignore_attr = TRUE)
})

test_that("Results from the down knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[3]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[9]]))
  expect_equal(r, reference_data[[9]], ignore_attr = TRUE)
})

# Testing Overexpression Signatures

## Testing the results for the any filtered signature
test_that("Results from the up knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[4]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[10]]))
  expect_equal(r, reference_data[[10]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[4]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[11]]))
  expect_equal(r, reference_data[[11]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[4]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[12]]))
  expect_equal(r, reference_data[[12]], ignore_attr = TRUE)
})

## Testing the results for the up filtered signature
test_that("Results from the up knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[5]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[13]]))
  expect_equal(r, reference_data[[13]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[5]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[14]]))
  expect_equal(r, reference_data[[14]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[5]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[15]]))
  expect_equal(r, reference_data[[15]], ignore_attr = TRUE)
})

## Testing the results for the down filtered signature
test_that("Results from the down knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[6]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[16]]))
  expect_equal(r, reference_data[[16]], ignore_attr = TRUE)
})

test_that("Results from the down knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[6]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[17]]))
  expect_equal(r, reference_data[[17]], ignore_attr = TRUE)
})

test_that("Results from the down knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[6]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[18]]))
  expect_equal(r, reference_data[[18]], ignore_attr = TRUE)
})

# Testing Chemical Perturbagen Signatures

## Testing the results for the any filtered signature
test_that("Results from the up knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[7]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[19]]))
  expect_equal(r, reference_data[[19]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[7]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[20]]))
  expect_equal(r, reference_data[[20]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[7]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[21]]))
  expect_equal(r, reference_data[[21]], ignore_attr = TRUE)
})


## Testing the results for the up filtered signature
test_that("Results from the up knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[8]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[22]]))
  expect_equal(r, reference_data[[22]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[8]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[23]]))
  expect_equal(r, reference_data[[23]], ignore_attr = TRUE)
})

test_that("Results from the up knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[8]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[24]]))
  expect_equal(r, reference_data[[24]], ignore_attr = TRUE)
})

## Testing the results for the down filtered signature
test_that("Results from the down knockdown signature are correct with respect to OE",{
  r <- get_concordants(input_data[[9]], "OE")
  expect_equal(nrow(r), nrow(reference_data[[25]]))
  expect_equal(r, reference_data[[25]], ignore_attr = TRUE)
})

test_that("Results from the down knockdown signature are correct with respect to KD",{
  r <- get_concordants(input_data[[9]], "KD")
  expect_equal(nrow(r), nrow(reference_data[[26]]))
  expect_equal(r, reference_data[[26]], ignore_attr = TRUE)
})

test_that("Results from the down knockdown signature are correct with respect to CP",{
  r <- get_concordants(input_data[[9]], "CP")
  expect_equal(nrow(r), nrow(reference_data[[27]]))
  expect_equal(r, reference_data[[27]], ignore_attr = TRUE)
})

