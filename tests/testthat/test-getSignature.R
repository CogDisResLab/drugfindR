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
vcr::use_cassette("test_valid", {
    test_that("getSignature has correct output for a valid signature", {
        validSignature <- getSignature("LINCSKD_28")
        expect_identical(nrow(validSignature), 978L)
        expect_named(validSignature, signatureColNames())
        expect_false(
            any(purrr::flatten_lgl(purrr::map(validSignature, is.na)))
        )
        expect_equal(validSignature, exampleSignature(), tolerance = 1e-12)
    })
})
