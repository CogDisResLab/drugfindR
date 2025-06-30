# Test the `getConcordants` function

# Test invalid signature

test_that("Function errors if it receives an error response", {
    webmockr::stub_request(
        "post", "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze"
    ) |>
        webmockr::to_return(status = 500L)
    webmockr::httr2_mock()
    expect_error(getConcordants(exampleSignature()))
    webmockr::httr2_mock(FALSE)
})


# Test valid signature

# test_that("getConcordants correct value", {
#     inputSignature <- exampleSignature() |>
#         filterSignature(threshold = 1.0)
#     concordantsList <- getConcordants(inputSignature, "CP")
#     expect_s3_class(concordantsList, "tbl_df")
#     expect_equal(concordantsList, concordantsCp(), tolerance = 1e-12)
#     expect_identical(ncol(concordantsList), 8L)
#     expect_identical(nrow(concordantsList), 14337L)
#     expect_identical(unique(concordantsList[["sig_direction"]]), "Any")
# })

test_that("getConcordants handles getting conrdants", {
    dfcd <- drugfindRCoreData.fixture() |>
        setFilterThresholds(1L) |>
        setOutputLibrary("CP") |>
        filterSignature() |>
        getConcordants()

    expect_s4_class(dfcd, "drugfindRResult")
})
