# Test the invalid signature

vcr::use_cassette(
    "test_invalid",
    {
        test_that("everything NA for invalid signature", {
            expect_error(getSignature("LINCS_INV"))
        })
    }
)

# Testing the retrieved signature
test_that("correct number of rows for the retrieved signature", {
    validSignature <- getSignature("LINCSKD_28")
    expect_identical(nrow(validSignature), 978L)
})

test_that("correct columns for the retrieved signature", {
    validSignature <- getSignature("LINCSKD_28")
    expect_named(validSignature, signatureColNames())
})

test_that("nothing NA for knockdown signature", {
    validSignature <- getSignature("LINCSKD_28")
    expect_false(
        any(purrr::flatten_lgl(purrr::map(validSignature, is.na)))
    )
})

test_that("correct content for the knockdown signature", {
    validSignature <- getSignature("LINCSKD_28")
    expect_equal(validSignature, exampleSignature(), tolerance = 1e-12)
})
