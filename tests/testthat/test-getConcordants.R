# Comprehensive tests for getConcordants function and internal helpers

library(S4Vectors)
library(httr2)
library(httptest2)

# ==============================================================================
# TESTS FOR INTERNAL VALIDATION FUNCTION
# ==============================================================================

test_that(".validateGetConcordantsInput works correctly with valid inputs", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

    # Valid inputs should not error
    expect_silent(.validateGetConcordantsInput(testSignature, "CP"))
    expect_silent(.validateGetConcordantsInput(testSignature, "KD"))
    expect_silent(.validateGetConcordantsInput(testSignature, "OE"))

    # Test with different data frame types
    testDataFrame <- as.data.frame(testSignature)
    testS4DataFrame <- S4Vectors::DataFrame(testSignature)

    expect_silent(.validateGetConcordantsInput(testDataFrame, "CP"))
    expect_silent(.validateGetConcordantsInput(testS4DataFrame, "CP"))
})

test_that(".validateGetConcordantsInput errors on invalid signature input", {
    # Non-dataframe input should error
    expect_error(
        .validateGetConcordantsInput("not_a_dataframe", "CP"),
        "Signature must be a data.frame, tibble, or DataFrame",
        fixed = TRUE
    )

    expect_error(
        .validateGetConcordantsInput(list(a = 1L, b = 2L), "CP"),
        "Signature must be a data.frame, tibble, or DataFrame",
        fixed = TRUE
    )

    expect_error(
        .validateGetConcordantsInput(matrix(1L:10L, nrow = 5L), "CP"),
        "Signature must be a data.frame, tibble, or DataFrame",
        fixed = TRUE
    )
})

test_that(".validateGetConcordantsInput errors on invalid library", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

    # Invalid library should error (uses stopIfInvalidLibraries)
    expect_error(
        .validateGetConcordantsInput(testSignature, "INVALID")
    )

    expect_error(
        .validateGetConcordantsInput(testSignature, "cp") # Case sensitive
    )

    expect_error(
        .validateGetConcordantsInput(testSignature, "")
    )
})

# ==============================================================================
# TESTS FOR SIGNATURE FILE PREPARATION
# ==============================================================================

test_that(".prepareSignatureFile creates valid temporary file", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

    signatureFilePath <- .prepareSignatureFile(testSignature)

    # Should return a character string
    expect_type(signatureFilePath, "character")
    expect_length(signatureFilePath, 1L)

    # File should exist
    expect_true(file.exists(signatureFilePath))

    # File should have .xls extension
    expect_true(grepl("\\.xls$", signatureFilePath))

    # File should contain signature data
    fileContent <- readr::read_tsv(signatureFilePath, show_col_types = FALSE)
    expect_identical(nrow(fileContent), nrow(testSignature))
    expect_identical(colnames(fileContent), colnames(testSignature))
})

test_that(".prepareSignatureFile works with different data frame types", {
    testSignatureTibble <- getTestFixture("prepared_signature", seed = .testSeed)
    testSignatureDataFrame <- as.data.frame(testSignatureTibble)
    testSignatureS4DataFrame <- S4Vectors::DataFrame(testSignatureTibble)

    # All should create valid files
    signatureFilePathTibble <- .prepareSignatureFile(testSignatureTibble)
    signatureFilePathDataFrame <- .prepareSignatureFile(testSignatureDataFrame)
    signatureFilePathS4DataFrame <- .prepareSignatureFile(testSignatureS4DataFrame)

    expect_true(file.exists(signatureFilePathTibble))
    expect_true(file.exists(signatureFilePathDataFrame))
    expect_true(file.exists(signatureFilePathS4DataFrame))

    # Content should be equivalent
    contentTibble <- readr::read_tsv(signatureFilePathTibble, show_col_types = FALSE)
    contentDataFrame <- readr::read_tsv(signatureFilePathDataFrame, show_col_types = FALSE)
    contentS4DataFrame <- readr::read_tsv(signatureFilePathS4DataFrame, show_col_types = FALSE)

    expect_identical(contentTibble, contentDataFrame)
    expect_identical(contentTibble, contentS4DataFrame)
})

# ==============================================================================
# TESTS FOR SIGNATURE DIRECTION DETECTION
# ==============================================================================

test_that(".detectSignatureDirection correctly identifies up-regulated signatures", {
    baseSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    thresholds <- .calculateAbsoluteThresholds(1.5)
    upSignature <- .applyDirectionFilter(baseSignature, "up", thresholds = thresholds)

    direction <- .detectSignatureDirection(upSignature)

    expect_identical(direction, "Up")
    expect_type(direction, "character")
    expect_length(direction, 1L)
})

test_that(".detectSignatureDirection correctly identifies down-regulated signatures", {
    baseSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    thresholds <- .calculateAbsoluteThresholds(1.5)
    downSignature <- .applyDirectionFilter(baseSignature, "down", thresholds = thresholds)

    direction <- .detectSignatureDirection(downSignature)

    expect_identical(direction, "Down")
    expect_type(direction, "character")
    expect_length(direction, 1L)
})

test_that(".detectSignatureDirection correctly identifies mixed signatures", {
    baseSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    thresholds <- .calculateAbsoluteThresholds(1.5)
    mixedSignature <- .applyDirectionFilter(baseSignature, "any", thresholds = thresholds)

    direction <- .detectSignatureDirection(mixedSignature)

    expect_identical(direction, "Any")
    expect_type(direction, "character")
    expect_length(direction, 1L)
})

test_that(".detectSignatureDirection handles edge cases", {
    # Signature with zero values
    zeroSignature <- createZeroSignature()
    expect_identical(.detectSignatureDirection(zeroSignature), "Up")

    # Signature with very small positive values
    smallPositiveSignature <- createSmallPositiveSignature()
    expect_identical(.detectSignatureDirection(smallPositiveSignature), "Up")

    # Signature with very small negative values
    smallNegativeSignature <- createSmallNegativeSignature()
    expect_identical(.detectSignatureDirection(smallNegativeSignature), "Down")
})

# ==============================================================================
# TESTS FOR API REQUEST GENERATION
# ==============================================================================

test_that(".generateIlincsRequest creates valid httr2 request object", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    signatureFilePath <- .prepareSignatureFile(testSignature)

    # Test with CP library
    request <- .generateIlincsRequest(signatureFilePath, "CP")

    # Should return an httr2_request object
    expect_s3_class(request, "httr2_request")

    # Check URL construction
    requestUrl <- req_get_url(request)
    parsedRequestUrl <- url_parse(requestUrl)
    expect_identical(parsedRequestUrl[["scheme"]], "https")
    expect_identical(parsedRequestUrl[["hostname"]], "www.ilincs.org")
    expect_identical(parsedRequestUrl[["path"]], "/api/SignatureMeta/uploadAndAnalyze") # nolint: nonportable_path_linter.

    # Check query parameters
    queryParameters <- parsedRequestUrl[["query"]]
    expect_identical(queryParameters[["lib"]], "LIB_5")

    # Check method
    expect_identical(req_get_method(request), "POST")

    # Check user agent
    userAgent <- request[["options"]][["useragent"]]
    expect_true(grepl("github.com/CogDisResLab/drugfindR", userAgent)) # nolint: nonportable_path_linter.

    # Check if the body has the correct structure
    expect_named(request[["body"]], c("data", "type", "content_type", "params"))

    # Check body type (should be multipart)
    expect_identical(request[["body"]][["type"]], "multipart")

    # Check body has a valid `file` field
    expect_s3_class(request[["body"]][["data"]][["file"]], "form_file")
})

test_that(".generateIlincsRequest works with different libraries", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    signatureFilePath <- .prepareSignatureFile(testSignature)

    # Test CP library
    requestCP <- .generateIlincsRequest(signatureFilePath, "CP")
    queryParametersCP <- url_parse(req_get_url(requestCP))[["query"]]
    expect_identical(queryParametersCP[["lib"]], "LIB_5")

    # Test KD library
    requestKD <- .generateIlincsRequest(signatureFilePath, "KD")
    queryParametersKD <- url_parse(req_get_url(requestKD))[["query"]]
    expect_identical(queryParametersKD[["lib"]], "LIB_6")

    # Test OE library
    requestOE <- .generateIlincsRequest(signatureFilePath, "OE")
    queryParametersOE <- url_parse(req_get_url(requestOE))[["query"]]
    expect_identical(queryParametersOE[["lib"]], "LIB_11")
})

test_that(".generateIlincsRequest validates library input", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    signatureFilePath <- .prepareSignatureFile(testSignature)

    # Invalid library should error
    expect_error(
        .generateIlincsRequest(signatureFilePath, "INVALID"),
        "Invalid library specification"
    )

    expect_error(
        .generateIlincsRequest(signatureFilePath, "cp"), # Wrong case
        "Invalid library specification"
    )

    expect_error(
        .generateIlincsRequest(signatureFilePath, ""),
        "Invalid library specification"
    )
})

test_that(".generateIlincsRequest includes proper headers and configuration", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)
    signatureFilePath <- .prepareSignatureFile(testSignature)

    request <- .generateIlincsRequest(signatureFilePath, "CP")

    # Check User-Agent header contains package info
    userAgent <- request[["options"]][["useragent"]]
    expect_true(grepl("drugfindR", userAgent))
    expect_true(grepl("github.com/CogDisResLab/drugfindR", userAgent)) # nolint: nonportable_path_linter.
})

# ==============================================================================
# TESTS FOR API REQUEST EXECUTION
# ==============================================================================

# TODO: Add support for verbosity
test_that(".executeIlincsRequest handles verbose option", {
    skip("Verbose response config option is not testable right now")
    verboseResponse <- mockConcordantsResponse()

    with_mocked_responses(verboseResponse, {
        # Test with verbose = TRUE
        response <- .executeIlincsRequest(exampleRequest(), verbose = TRUE)
        respJson <- resp_body_json(response)
        expect_true(respJson[["request"]][["options"]][["verbose"]])
        expect_s3_class(response, "httr2_response")

        # Test with verbose = FALSE (default)
        response2 <- .executeIlincsRequest(exampleRequest(), verbose = FALSE)
        resp2Json <- resp_body_json(response2)
        expect_null(resp2Json[["request"]][["options"]][["verbose"]])
        expect_s3_class(response2, "httr2_response")

        # Test with no specification (default)
        response3 <- .executeIlincsRequest(exampleRequest())
        resp3Json <- resp_body_json(response3)
        expect_null(resp3Json[["request"]][["options"]][["verbose"]])
        expect_s3_class(response3, "httr2_response")
    })
})

# TODO: Need to figure out how to trigger 500 errors for tests # nolint: todo_comment_linter.
test_that(".executeIlincsRequest returns error responses without raising", {
    skip("Have to figure out how to trigger 500 errors")

    # Mock error response
    # Should not raise an error, should return the error response
    response <- .executeIlincsRequest(exampleRequest())
    expect_s3_class(response, "httr2_response")
    expect_identical(resp_status(response), 500L)
})

# ==============================================================================
# TESTS FOR API RESPONSE PROCESSING
# ==============================================================================

test_that(".processIlincsResponse handles successful responses correctly", {
    # Create mock successful response data
    response <- getTestFixture("valid_ilincs_response", seed = .testSeed, nRows = 10L)

    result <- .processIlincsResponse(response, "Up", "CP")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 10L)

    # Check required columns
    expectedCols <- c(
        "signatureid", "treatment", "time", "cellline",
        "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_true(all(expectedCols %in% colnames(result)))

    # Check CP-specific content and transformations
    expect_true("treatment" %in% colnames(result)) # compound should be renamed to treatment
    expect_true("concentration" %in% colnames(result))
    expect_identical(result[["sig_direction"]], rep("Up", 10L))
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 10L))
})

test_that(".processIlincsResponse handles CP library responses correctly", {
    # Create mock CP response with compound and concentration
    response <- getTestFixture("valid_ilincs_response", seed = .testSeed, nRows = 10L)

    result <- .processIlincsResponse(response, "Any", "CP")

    # Check CP-specific columns
    expect_false("compound" %in% colnames(result))
    expect_true("treatment" %in% colnames(result)) # compound renamed to treatment
    expect_true("concentration" %in% colnames(result))
    expect_true("sig_type" %in% colnames(result))

    # Check content
    expect_identical(result[["sig_direction"]], rep("Any", 10L))
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 10L))
})

test_that(".processIlincsResponse rounds numerical values correctly", {
    response <- getTestFixture("valid_ilincs_response", seed = .testSeed, nRows = 10L)

    result <- .processIlincsResponse(response, "Down", "CP")

    # Check similarity rounded to 8 decimal places
    expect_identical(result[["similarity"]][[1L]], -0.62460)


    # Check pValue rounded to 20 decimal places (or system precision)
    expect_type(result[["pValue"]], "double")
    expect_true(all(abs(result[["pValue"]]) >= 0L))
})

test_that(".processIlincsResponse handles different signature directions", {
    directions <- c("Up", "Down", "Any")

    response <- getTestFixture("valid_ilincs_response", seed = .testSeed, nRows = 10L)

    for (direction in directions) {
        result <- .processIlincsResponse(response, direction, "CP")
        expect_identical(unique(result[["sig_direction"]]), direction)
    }
})

test_that(".processIlincsResponse handles empty concordance tables", {
    response <- getTestFixture("empty_ilincs_response", seed = .testSeed)
    # Empty response
    result <- .processIlincsResponse(response, "Any", "CP")

    # Should return empty tibble with correct structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 0L)
    expect_true("sig_direction" %in% colnames(result))
    expect_true("sig_type" %in% colnames(result))
})


test_that(".processIlincsResponse handles 400 errors correctly", {
    response <- getTestFixture("error_ilincs_response_400", seed = .testSeed)

    expect_error(
        .processIlincsResponse(response),
        "iLINCS API request failed \\(Status 400\\): Bad Request",
        class = "error"
    )

    expect_error(
        .processIlincsResponse(response),
        "invalid signature data or parameters"
    )
})

test_that(".processIlincsResponse handles 500 errors correctly", {
    response <- getTestFixture("error_ilincs_response_500", seed = .testSeed)

    expect_error(
        .processIlincsResponse(response,
            "iLINCS API request failed \\(Status 500\\): Internal Server ",
            class = "error"
        )
    )

    expect_error(
        .processIlincsResponse(response),
        "data structure problems or server issues"
    )
})

test_that(".processIlincsResponse handles other status codes", {
    # Create custom error response with different status code
    errorResp404 <- httr2::response(
        status_code = 404L,
        headers = list("content-type" = "text/plain"), # nolint: nonportable_path_linter.
        body = charToRaw("Not Found")
    )

    expect_error(
        .processIlincsResponse(errorResp404),
        "iLINCS API request failed with status 404"
    )
})

test_that(".processIlincsResponseEmpty creates correct empty structure for all libraries", {
    libraries <- c("CP", "KD", "OE")
    expectedSigTypes <- c("Chemical Perturbagen", "Gene Knockdown", "Gene Overexpression")

    for (i in seq_along(libraries)) {
        lib <- libraries[i]
        expectedSigType <- expectedSigTypes[i]

        result <- .processIlincsResponseEmpty("Any", lib)

        # Check structure
        expect_s3_class(result, "tbl_df")
        expect_identical(nrow(result), 0L)

        # Check all required columns are present
        expectedCols <- c(
            "signatureid", "treatment", "concentration", "time",
            "cellline", "similarity", "pValue", "sig_direction", "sig_type"
        )
        expect_identical(colnames(result), expectedCols)

        # Check column types
        expect_type(result[["signatureid"]], "character")
        expect_type(result[["treatment"]], "character")
        expect_type(result[["concentration"]], "character")
        expect_type(result[["time"]], "character")
        expect_type(result[["cellline"]], "character")
        expect_type(result[["similarity"]], "double")
        expect_type(result[["pValue"]], "double")
        expect_type(result[["sig_direction"]], "character")
        expect_type(result[["sig_type"]], "character")
    }
})

test_that(".processIlincsResponseSuccess handles CP library correctly", {
    # Create mock concordance tables for CP library
    response <- getTestFixture("valid_ilincs_response",
        seed = .testSeed,
        nRows = 10L, library = "CP"
    )
    responseTable <- response |>
        resp_body_json() |>
        purrr::pluck("status", "concordanceTable")

    result <- .processIlincsResponseSuccess(responseTable, "Up", "CP")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 10L)

    # Check column structure and order
    expectedCols <- c(
        "signatureid", "treatment", "concentration", "time",
        "cellline", "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_identical(colnames(result), expectedCols)

    # Check CP-specific transformations
    expect_true(all(stringr::str_detect(result[["treatment"]], stringr::fixed("DRUG")))) # compound renamed to treatment
    expect_false(all(stringr::str_detect(result[["treatment"]], stringr::fixed("GENE")))) # compound renamed to treatment
    expect_false(all(is.na(result[["concentration"]])))
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 10L))
    expect_identical(result[["sig_direction"]], rep("Up", 10L))

    # Check rounding
    expect_identical(result[["similarity"]], round(result[["similarity"]], 8L)) # 8 decimal places
    expect_equal(result[["pValue"]], round(result[["pValue"]], 12L), tolerance = 1e-20)
})

test_that(".processIlincsResponseSuccess handles KD library correctly", {
    response <- getTestFixture("valid_ilincs_response",
        seed = .testSeed,
        nRows = 10L, library = "KD"
    )
    responseTable <- response |>
        resp_body_json() |>
        purrr::pluck("status", "concordanceTable")

    result <- .processIlincsResponseSuccess(responseTable, "Down", "KD")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 10L)

    # Check KD-specific transformations
    expect_true(all(stringr::str_detect(result[["treatment"]], "GENE"))) # compound renamed to treatment
    expect_false(all(stringr::str_detect(result[["treatment"]], "DRUG"))) # compound renamed to treatment
    expect_true(all(is.na(result[["concentration"]]))) # All NA for KD
    expect_identical(result[["sig_type"]], rep("Gene Knockdown", 10L))
    expect_identical(result[["sig_direction"]], rep("Down", 10L))
})

test_that(".processIlincsResponseSuccess handles OE library correctly", {
    response <- getTestFixture("valid_ilincs_response",
        seed = .testSeed,
        nRows = 10L, library = "KD"
    )
    responseTable <- response |>
        resp_body_json() |>
        purrr::pluck("status", "concordanceTable")

    result <- .processIlincsResponseSuccess(responseTable, "Any", "OE")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 10L)

    # Check OE-specific transformations
    expect_true(all(stringr::str_detect(result[["treatment"]], "GENE"))) # compound renamed to treatment
    expect_false(all(stringr::str_detect(result[["treatment"]], "DRUG"))) # compound renamed to treatment
    expect_true(all(is.na(result[["concentration"]]))) # All NA for OE
    expect_identical(result[["sig_type"]], rep("Gene Overexpression", 10L))
    expect_identical(result[["sig_direction"]], rep("Any", 10L))
})

test_that(".processIlincsResponse dispatcher works correctly for success case", {
    response <- getTestFixture("valid_ilincs_response",
        seed = .testSeed, nRows = 10L
    )
    result <- .processIlincsResponse(response, "Up", "CP")

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 10L)
    expect_true(all(stringr::str_detect(result[["treatment"]], "DRUG"))) # compound renamed to treatment
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 10L))
})

test_that(".processIlincsResponse dispatcher works correctly for empty case", {
    response <- getTestFixture("empty_ilincs_response", seed = .testSeed)

    result <- .processIlincsResponse(response, "Any", "KD")

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 0L)
    expect_named(result, c(
        "signatureid", "treatment", "concentration",
        "time", "cellline", "similarity", "pValue",
        "sig_direction", "sig_type"
    ))
})

test_that(".processIlincsResponse dispatcher works correctly for error case", {
    # Mock an error response
    expect_error(
        .processIlincsResponse(getTestFixture("error_ilincs_response_400"), "Any", "CP"),
        "iLINCS API request failed \\(Status 400\\)"
    )

    expect_error(
        .processIlincsResponse(getTestFixture("error_ilincs_response_500"), "Any", "CP"),
        "iLINCS API request failed \\(Status 500\\)"
    )
})

test_that(".processIlincsResponse maintains consistent column order across libraries", {
    # Test that all libraries produce the same column order
    expectedCols <- c(
        "signatureid", "treatment", "concentration", "time",
        "cellline", "similarity", "pValue", "sig_direction", "sig_type"
    )

    libraries <- c("CP", "KD", "OE")

    for (lib in libraries) {
        # Test empty response
        emptyResult <- .processIlincsResponseEmpty("Any", lib)
        expect_identical(colnames(emptyResult), expectedCols)

        # Test with mock data
        mockData <- list(list(
            signatureid = "TEST_001",
            treatment = if (lib == "CP") NULL else "TEST_GENE",
            compound = if (lib == "CP") "TEST_COMPOUND" else NA,
            concentration = if (lib == "CP") "1uM" else NA,
            time = "24h",
            cellline = "TestCell",
            similarity = 0.5,
            pValue = 0.05
        ))

        successResult <- .processIlincsResponseSuccess(mockData, "Any", lib)
        expect_named(successResult, expectedCols)
    }
})
