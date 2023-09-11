# Test the invalid signature
test_that("everything NA for invalid signature", {
    invalid_signature <- get_signature("LINCS_INV")

    expect_true(
        all(purrr::flatten_lgl(purrr::map(invalid_signature, is.na)))
    )
})

test_that("correct columns for the invalid signature", {
    invalid_signature <- get_signature("LINCS_INV")
    expect_named(invalid_signature, signature_col_names())
})

test_that("correct content for the invalid signature", {
    invalid_signature <- get_signature("LINCS_INV")
    expect_identical(invalid_signature, empty_signature())
})

# Testing the retrieved signature
test_that("correct number of rows for the retrieved signature", {
    valid_signature <- get_signature("LINCSKD_28")
    expect_identical(nrow(valid_signature), 978L)
})

test_that("correct columns for the retrieved signature", {
    valid_signature <- get_signature("LINCSKD_28")
    expect_named(valid_signature, signature_col_names())
})

test_that("nothing NA for knockdown signature", {
    valid_signature <- get_signature("LINCSKD_28")
    expect_false(
        any(purrr::flatten_lgl(purrr::map(valid_signature, is.na)))
    )
})

test_that("correct content for the knockdown signature", {
    valid_signature <- get_signature("LINCSKD_28")
    expect_equal(valid_signature, example_signature(), tolerance = 1e-12)
})
