# Test the `get_concordants` function

# Test Invalid inputs

test_that("Input signature must be a data frame or data frame like object", {
    expect_error(get_concordants("LINCSKD_28"))
})

test_that("Library must be one of 'OE', 'KD' or 'CP'", {
    expect_error(get_concordants(example_signature(), "INVALID"))
})

# Test invalid signature

test_that("Function errors if it receives an error response", {
    webmockr::stub_request("post", "http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze") %>%
        webmockr::to_return(status = 500L)
    webmockr::httr_mock()
    expect_error(get_concordants(example_signature()))
    webmockr::httr_mock(FALSE)
})


# Test valid signature

test_that("get_concordants correct value", {
    input_signature <- example_signature() |>
        filter_signature(threshold = 1.0)
    concordants_list <- get_concordants(input_signature, "CP", "any")
    expect_s3_class(concordants_list, "tbl_df")
    expect_equal(concordants_list, concordants_cp(), tolerance = 1e-12)
    expect_identical(ncol(concordants_list), 8L)
    expect_identical(nrow(concordants_list), 14337L)
    expect_identical(unique(concordants_list[["sig_direction"]]), "any")
})
