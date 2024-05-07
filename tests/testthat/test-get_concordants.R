# Test the `getConcordants` function

# Test Invalid inputs

test_that("Input signature must be a data frame or data frame like object", {
    expect_error(getConcordants("LINCSKD_28"))
})

test_that("Library must be one of 'OE', 'KD' or 'CP'", {
    expect_error(getConcordants(exampleSignature(), "INVALID"))
})

# Test invalid signature

test_that("Function errors if it receives an error response", {
    webmockr::stub_request(
        "post", "http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze"
    ) |>
        webmockr::to_return(status = 500L)
    webmockr::httr_mock()
    expect_error(getConcordants(exampleSignature()))
    webmockr::httr_mock(FALSE)
})


# Test valid signature

test_that("getConcordants correct value", {
    inputSignature <- exampleSignature() |>
        filterSignature(threshold = 1.0)
    concordantsList <- getConcordants(inputSignature, "CP", "any")
    expect_s3_class(concordantsList, "tbl_df")
    expect_equal(concordantsList, concordantsCp(), tolerance = 1e-12)
    expect_identical(ncol(concordantsList), 8L)
    expect_identical(nrow(concordantsList), 14337L)
    expect_identical(unique(concordantsList[["sig_direction"]]), "any")
})
