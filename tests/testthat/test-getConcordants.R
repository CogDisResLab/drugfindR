# Comprehensive tests for getConcordants function and internal helpers

library(tibble)
library(S4Vectors)
library(vcr)
library(httr2)

# ==============================================================================
# TESTS FOR INTERNAL VALIDATION FUNCTION
# ==============================================================================

test_that(".validateGetConcordantsInput works correctly with valid inputs", {
    testSig <- exampleSignature()

    # Valid inputs should not error
    expect_silent(.validateGetConcordantsInput(testSig, "CP"))
    expect_silent(.validateGetConcordantsInput(testSig, "KD"))
    expect_silent(.validateGetConcordantsInput(testSig, "OE"))

    # Test with different data frame types
    testDf <- as.data.frame(testSig)
    testDataFrame <- S4Vectors::DataFrame(testSig)

    expect_silent(.validateGetConcordantsInput(testDf, "CP"))
    expect_silent(.validateGetConcordantsInput(testDataFrame, "CP"))
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
    testSig <- exampleSignature()

    # Invalid library should error (uses stopIfInvalidLibraries)
    expect_error(
        .validateGetConcordantsInput(testSig, "INVALID")
    )

    expect_error(
        .validateGetConcordantsInput(testSig, "cp") # Case sensitive
    )

    expect_error(
        .validateGetConcordantsInput(testSig, "")
    )
})

# ==============================================================================
# TESTS FOR SIGNATURE FILE PREPARATION
# ==============================================================================

test_that(".prepareSignatureFile creates valid temporary file", {
    testSig <- exampleSignature()

    filePath <- .prepareSignatureFile(testSig)

    # Should return a character string
    expect_type(filePath, "character")
    expect_length(filePath, 1L)

    # File should exist
    expect_true(file.exists(filePath))

    # File should have .xls extension
    expect_true(grepl("\\.xls$", filePath))

    # File should contain signature data
    fileContent <- readr::read_tsv(filePath, show_col_types = FALSE)
    expect_identical(nrow(fileContent), nrow(testSig))
    expect_identical(colnames(fileContent), colnames(testSig))
})

test_that(".prepareSignatureFile works with different data frame types", {
    testDataTibble <- exampleSignature()
    testDataDf <- as.data.frame(testDataTibble)
    testDataDataFrame <- S4Vectors::DataFrame(testDataTibble)

    # All should create valid files
    filePathTibble <- .prepareSignatureFile(testDataTibble)
    filePathDf <- .prepareSignatureFile(testDataDf)
    filePathDataFrame <- .prepareSignatureFile(testDataDataFrame)

    expect_true(file.exists(filePathTibble))
    expect_true(file.exists(filePathDf))
    expect_true(file.exists(filePathDataFrame))

    # Content should be equivalent
    contentTibble <- readr::read_tsv(filePathTibble, show_col_types = FALSE)
    contentDf <- readr::read_tsv(filePathDf, show_col_types = FALSE)
    contentDataFrame <- readr::read_tsv(filePathDataFrame, show_col_types = FALSE)

    expect_identical(contentTibble, contentDf)
    expect_identical(contentTibble, contentDataFrame)
})

# ==============================================================================
# TESTS FOR SIGNATURE DIRECTION DETECTION
# ==============================================================================

test_that(".detectSignatureDirection correctly identifies up-regulated signatures", {
    upSig <- exampleSignatureUpFilter()

    direction <- .detectSignatureDirection(upSig)

    expect_identical(direction, "Up")
    expect_type(direction, "character")
    expect_length(direction, 1L)
})

test_that(".detectSignatureDirection correctly identifies down-regulated signatures", {
    downSig <- exampleSignatureDownFilter()

    direction <- .detectSignatureDirection(downSig)

    expect_identical(direction, "Down")
    expect_type(direction, "character")
    expect_length(direction, 1L)
})

test_that(".detectSignatureDirection correctly identifies mixed signatures", {
    mixedSig <- exampleSignatureAnyFilter() # Has both positive and negative values

    direction <- .detectSignatureDirection(mixedSig)

    expect_identical(direction, "Any")
    expect_type(direction, "character")
    expect_length(direction, 1L)
})

test_that(".detectSignatureDirection handles edge cases", {
    # Signature with zero values
    zeroSig <- tibble::tibble(
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = 0.0
    )
    expect_identical(.detectSignatureDirection(zeroSig), "Up")

    # Signature with very small positive values
    smallPosSig <- tibble::tibble(
        Name_GeneSymbol = c("GENE1", "GENE2"),
        Value_LogDiffExp = c(0.001, 0.002)
    )
    expect_identical(.detectSignatureDirection(smallPosSig), "Up")

    # Signature with very small negative values
    smallNegSig <- tibble::tibble(
        Name_GeneSymbol = c("GENE1", "GENE2"),
        Value_LogDiffExp = c(-0.001, -0.002)
    )
    expect_identical(.detectSignatureDirection(smallNegSig), "Down")
})

# ==============================================================================
# TESTS FOR API REQUEST GENERATION
# ==============================================================================

test_that(".generateIlincsRequest creates valid httr2 request object", {
    testSig <- exampleSignature()
    filePath <- .prepareSignatureFile(testSig)

    # Test with CP library
    request <- .generateIlincsRequest(filePath, "CP")

    # Should return an httr2_request object
    expect_s3_class(request, "httr2_request")

    requestUrl <- req_get_url(request)

    # Check URL construction
    expect_true(grepl("ilincs.org/api", requestUrl))
    expect_true(grepl("SignatureMeta/uploadAndAnalyze", requestUrl))

    # Check query parameters
    expect_true(grepl("lib=LIB_5", requestUrl))

    # Check method
    expect_identical(req_get_method(request), "POST")

    # Check user agent
    expect_true(grepl("drugfindR", request[["options"]][["useragent"]]))

    # Check if the body has the correct structure
    expect_named(request[["body"]], c("data", "type", "content_type", "params"))

    # Check body type (should be multipart)
    expect_identical(request[["body"]][["type"]], "multipart")

    # Check body has a valid `file` field
    expect_s3_class(request[["body"]][["data"]][["file"]], "form_file")
})

test_that(".generateIlincsRequest works with different libraries", {
    testSig <- exampleSignature()
    filePath <- .prepareSignatureFile(testSig)

    # Test CP library
    requestCP <- .generateIlincsRequest(filePath, "CP")
    expect_true(grepl("lib=LIB_5", requestCP[["url"]]))

    # Test KD library
    requestKD <- .generateIlincsRequest(filePath, "KD")
    expect_true(grepl("lib=LIB_6", requestKD[["url"]]))

    # Test OE library
    requestOE <- .generateIlincsRequest(filePath, "OE")
    expect_true(grepl("lib=LIB_11", requestOE[["url"]]))
})

test_that(".generateIlincsRequest validates library input", {
    testSig <- exampleSignature()
    filePath <- .prepareSignatureFile(testSig)

    # Invalid library should error
    expect_error(
        .generateIlincsRequest(filePath, "INVALID"),
        "Invalid library specification"
    )

    expect_error(
        .generateIlincsRequest(filePath, "cp"), # Wrong case
        "Invalid library specification"
    )

    expect_error(
        .generateIlincsRequest(filePath, ""),
        "Invalid library specification"
    )
})

test_that(".generateIlincsRequest includes proper headers and configuration", {
    testSig <- exampleSignature()
    filePath <- .prepareSignatureFile(testSig)

    request <- .generateIlincsRequest(filePath, "CP")

    # Check User-Agent header contains package info
    userAgent <- request[["options"]][["useragent"]]
    expect_true(grepl("drugfindR", userAgent))
    expect_true(grepl("github.com/CogDisResLab/drugfindR", userAgent))

    # Check URL structure
    expect_true(grepl("^https://", request[["url"]]))
    expect_true(grepl("ilincs\\.org", request[["url"]]))
})

# ==============================================================================
# TESTS FOR API REQUEST EXECUTION
# ==============================================================================

test_that(".executeIlincsRequest configures request properly", {
    successResponse <- mockResponse(list(item = "Nothing"), 200L)

    with_mocked_responses(successResponse, {
        resp <- .executeIlincsRequest(exampleRequest())
    })

    expect_s3_class(resp, "httr2_response")
})

test_that(".executeIlincsRequest handles verbose option", {
    verboseResponse <- mockResponse(list(), 200L)

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

test_that(".executeIlincsRequest returns error responses without raising", {
    errorResponse <- mockResponse(list(), 500L)

    # Mock error response
    with_mocked_responses(errorResponse, {
        # Should not raise an error, should return the error response
        response <- .executeIlincsRequest(exampleRequest())
        expect_s3_class(response, "httr2_response")
        expect_identical(resp_status(response), 500L)
    })
})

# ==============================================================================
# TESTS FOR API RESPONSE PROCESSING
# ==============================================================================

test_that(".processIlincsResponse handles successful responses correctly", {
    # Create mock successful response data

    result <- .processIlincsResponse(exampleResponse(), "Up", "CP")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 6L)

    # Check required columns
    expectedCols <- c(
        "signatureid", "treatment", "time", "cellline",
        "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_true(all(expectedCols %in% colnames(result)))

    # Check CP-specific content and transformations
    expect_true("treatment" %in% colnames(result)) # compound should be renamed to treatment
    expect_true("concentration" %in% colnames(result))
    expect_identical(result[["sig_direction"]], rep("Up", 6L))
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 6L))
})

test_that(".processIlincsResponse handles CP library responses correctly", {
    # Create mock CP response with compound and concentration

    result <- .processIlincsResponse(exampleResponse(), "Any", "CP")

    # Check CP-specific columns
    expect_true("treatment" %in% colnames(result)) # compound renamed to treatment
    expect_true("concentration" %in% colnames(result))
    expect_true("sig_type" %in% colnames(result))

    # Check content
    expect_identical(result[["sig_direction"]], rep("Any", 6L))
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 6L))
})

test_that(".processIlincsResponse rounds numerical values correctly", {
    result <- .processIlincsResponse(exampleResponse(), "Down", "CP")

    # Check similarity rounded to 8 decimal places
    expect_identical(result[["similarity"]][[1L]], 0.54640)
    expect_true(all(abs(result[["similarity"]]) >= 0.2))


    # Check pValue rounded to 20 decimal places (or system precision)
    expect_type(result[["pValue"]], "double")
    expect_true(all(abs(result[["pValue"]]) >= 0L))
})

test_that(".processIlincsResponse handles different signature directions", {
    directions <- c("Up", "Down", "Any")

    for (direction in directions) {
        result <- .processIlincsResponse(exampleResponse(), direction, "CP")
        expect_identical(unique(result[["sig_direction"]]), direction)
    }
})

test_that(".processIlincsResponse handles empty concordance tables", {
    # Empty response
    result <- .processIlincsResponse(emptyResponse(), "Any", "CP")

    # Should return empty tibble with correct structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 0L)
    expect_true("sig_direction" %in% colnames(result))
    expect_true("sig_type" %in% colnames(result))
})

# ==============================================================================
# TESTS FOR NEW PROCESSILINCSRESPONSE FAMILY OF FUNCTIONS
# ==============================================================================

test_that(".processIlincsResponseError handles 400 errors correctly", {
    errorResp400 <- errorResponse400()

    expect_error(
        .processIlincsResponseError(errorResp400),
        "iLINCS API request failed \\(Status 400\\): Bad Request",
        class = "error"
    )

    expect_error(
        .processIlincsResponseError(errorResp400),
        "invalid signature data or parameters"
    )
})

test_that(".processIlincsResponseError handles 500 errors correctly", {
    errorResp500 <- errorResponse500()

    expect_error(
        .processIlincsResponseError(errorResp500),
        "iLINCS API request failed \\(Status 500\\): Internal Server Error",
        class = "error"
    )

    expect_error(
        .processIlincsResponseError(errorResp500),
        "data structure problems or server issues"
    )
})

test_that(".processIlincsResponseError handles other status codes", {
    # Create custom error response with different status code
    errorResp404 <- httr2::response(
        status_code = 404L,
        headers = list("content-type" = "text/plain"),
        body = charToRaw("Not Found")
    )

    expect_error(
        .processIlincsResponseError(errorResp404),
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
    mockCpTables <- list(
        list(
            signatureid = "CPC_001",
            compound = "Aspirin",
            concentration = "10uM",
            time = "24h",
            cellline = "HeLa",
            similarity = 0.123456789,
            pValue = 0.001234567890123456789
        ),
        list(
            signatureid = "CPC_002",
            compound = "Ibuprofen",
            concentration = "5uM",
            time = "48h",
            cellline = "MCF7",
            similarity = -0.987654321,
            pValue = 0.045678901234567890123
        )
    )

    result <- .processIlincsResponseSuccess(mockCpTables, "Up", "CP")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 2L)

    # Check column structure and order
    expectedCols <- c(
        "signatureid", "treatment", "concentration", "time",
        "cellline", "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_identical(colnames(result), expectedCols)

    # Check CP-specific transformations
    expect_identical(result[["treatment"]], c("Aspirin", "Ibuprofen")) # compound renamed to treatment
    expect_identical(result[["concentration"]], c("10uM", "5uM"))
    expect_identical(result[["sig_type"]], c("Chemical Perturbagen", "Chemical Perturbagen"))
    expect_identical(result[["sig_direction"]], c("Up", "Up"))

    # Check rounding
    expect_identical(result[["similarity"]], c(0.12345679, -0.98765432)) # 8 decimal places
    expect_equal(result[["pValue"]], c(0.001234567890123456789, 0.045678901234567890123), tolerance = 1e-20)
})

test_that(".processIlincsResponseSuccess handles KD library correctly", {
    # Create mock concordance tables for KD library
    mockKdTables <- list(
        list(
            signatureid = "LINCSKD_001",
            treatment = "GENE1",
            time = "96h",
            cellline = "A549",
            similarity = 0.567890123,
            pValue = 0.012345678901234567890
        ),
        list(
            signatureid = "LINCSKD_002",
            treatment = "GENE2",
            time = "120h",
            cellline = "HepG2",
            similarity = -0.345678901,
            pValue = 0.034567890123456789012
        )
    )

    result <- .processIlincsResponseSuccess(mockKdTables, "Down", "KD")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 2L)

    # Check KD-specific transformations
    expect_identical(result[["treatment"]], c("GENE1", "GENE2"))
    expect_true(all(is.na(result[["concentration"]]))) # All NA for KD
    expect_identical(result[["sig_type"]], c("Gene Knockdown", "Gene Knockdown"))
    expect_identical(result[["sig_direction"]], c("Down", "Down"))
})

test_that(".processIlincsResponseSuccess handles OE library correctly", {
    # Create mock concordance tables for OE library
    mockOeTables <- list(
        list(
            signatureid = "LINCSOE_001",
            treatment = "GENE_X",
            time = "72h",
            cellline = "U2OS",
            similarity = 0.789012345,
            pValue = 0.023456789012345678901
        )
    )

    result <- .processIlincsResponseSuccess(mockOeTables, "Any", "OE")

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 1L)

    # Check OE-specific transformations
    expect_identical(result[["treatment"]], "GENE_X")
    expect_true(all(is.na(result[["concentration"]]))) # All NA for OE
    expect_identical(result[["sig_type"]], "Gene Overexpression")
    expect_identical(result[["sig_direction"]], "Any")
})

test_that(".processIlincsResponse dispatcher works correctly for success case", {
    # Mock a successful response with CP data
    result <- .processIlincsResponse(exampleResponse(), "Up", "CP")

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 6L)
    expect_identical(result[["treatment"]], c(
        "CHEMBL1271119", "MLS002264408",
        "VU0413807-2", "Flunisolide",
        "MLS002264499", "SMR000031368"
    ))
    expect_identical(result[["sig_type"]], rep("Chemical Perturbagen", 6L))
})

test_that(".processIlincsResponse dispatcher works correctly for empty case", {
    # Mock a response with empty concordance table (single "NA")

    result <- .processIlincsResponse(emptyResponse(), "Any", "KD")

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 0L)
    expect_identical(colnames(result), c(
        "signatureid", "treatment", "concentration",
        "time", "cellline", "similarity", "pValue",
        "sig_direction", "sig_type"
    ))
})

test_that(".processIlincsResponse dispatcher works correctly for error case", {
    # Mock an error response
    expect_error(
        .processIlincsResponse(errorResponse400(), "Any", "CP"),
        "iLINCS API request failed \\(Status 400\\)"
    )

    expect_error(
        .processIlincsResponse(errorResponse500(), "Any", "CP"),
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
            compound = if (lib == "CP") "TEST_COMPOUND" else NULL,
            concentration = if (lib == "CP") "1uM" else NULL,
            time = "24h",
            cellline = "TestCell",
            similarity = 0.5,
            pValue = 0.05
        ))

        successResult <- .processIlincsResponseSuccess(mockData, "Any", lib)
        expect_identical(colnames(successResult), expectedCols)
    }
})

# ==============================================================================
# VCR INTEGRATION TESTS - API-LEVEL CHECKS
# ==============================================================================

test_that("getConcordants works end-to-end with CP library via VCR", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()

    vcr::local_cassette("getConcordants_CP_integration")

    # Use a small signature for faster API calls
    testSig <- exampleSignature()[1L:5L, ]

    result <- getConcordants(testSig, ilincsLibrary = "CP")

    # Check result structure
    expect_s3_class(result, "tbl_df")
    expectedCols <- c(
        "signatureid", "treatment", "concentration", "time",
        "cellline", "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_identical(colnames(result), expectedCols)

    # Check CP-specific content
    if (nrow(result) > 0L) {
        expect_identical(unique(result[["sig_type"]]), "Chemical Perturbagen")
        expect_true(all(result[["sig_direction"]] %in% c("Up", "Down", "Any")))
        expect_true(all(is.numeric(result[["similarity"]])))
        expect_true(all(is.numeric(result[["pValue"]])))
    }
})

test_that("getConcordants works end-to-end with KD library via VCR", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()

    vcr::local_cassette("getConcordants_KD_integration")

    # Use a small signature for faster API calls
    testSig <- exampleSignature()[1L:5L, ]

    result <- getConcordants(testSig, ilincsLibrary = "KD")

    # Check result structure
    expect_s3_class(result, "tbl_df")
    expectedCols <- c(
        "signatureid", "treatment", "concentration", "time",
        "cellline", "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_identical(colnames(result), expectedCols)

    # Check KD-specific content
    if (nrow(result) > 0L) {
        expect_identical(unique(result[["sig_type"]]), "Gene Knockdown")
        expect_true(all(is.na(result[["concentration"]]))) # KD has no concentration
        expect_true(all(result[["sig_direction"]] %in% c("Up", "Down", "Any")))
        expect_true(all(is.numeric(result[["similarity"]])))
        expect_true(all(is.numeric(result[["pValue"]])))
    }
})

test_that("getConcordants works end-to-end with OE library via VCR", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()

    vcr::local_cassette("getConcordants_OE_integration")

    # Use a small signature for faster API calls
    testSig <- exampleSignature()[1L:5L, ]

    result <- getConcordants(testSig, ilincsLibrary = "OE")

    # Check result structure
    expect_s3_class(result, "tbl_df")
    expectedCols <- c(
        "signatureid", "treatment", "concentration", "time",
        "cellline", "similarity", "pValue", "sig_direction", "sig_type"
    )
    expect_identical(colnames(result), expectedCols)

    # Check OE-specific content
    if (nrow(result) > 0L) {
        expect_identical(unique(result[["sig_type"]]), "Gene Overexpression")
        expect_true(all(is.na(result[["concentration"]]))) # OE has no concentration
        expect_true(all(result[["sig_direction"]] %in% c("Up", "Down", "Any")))
        expect_true(all(is.numeric(result[["similarity"]])))
        expect_true(all(is.numeric(result[["pValue"]])))
    }
})

test_that("getConcordants handles different signature directions via VCR", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()

    vcr::local_cassette("getConcordants_directions_integration")

    # Test with up-regulated signature
    upSig <- exampleSignatureUpFilter()[1L:3L, ]
    upResult <- getConcordants(upSig, ilincsLibrary = "CP")

    if (nrow(upResult) > 0L) {
        expect_identical(unique(upResult[["sig_direction"]]), "Up")
    }

    # Test with down-regulated signature
    downSig <- exampleSignatureDownFilter()[1L:3L, ]
    downResult <- getConcordants(downSig, ilincsLibrary = "CP")

    if (nrow(downResult) > 0L) {
        expect_identical(unique(downResult[["sig_direction"]]), "Down")
    }

    # Test with mixed signature
    mixedSig <- exampleSignatureAnyFilter()[1L:3L, ]
    mixedResult <- getConcordants(mixedSig, ilincsLibrary = "CP")

    if (nrow(mixedResult) > 0L) {
        expect_identical(unique(mixedResult[["sig_direction"]]), "Any")
    }
})

test_that("getConcordants maintains input/output type consistency via VCR", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()

    vcr::local_cassette("getConcordants_types_integration")

    # Test with different input types
    testSig <- exampleSignature()[1L:3L, ]

    # Test with tibble input
    tibbleResult <- getConcordants(testSig, "CP")
    expect_s3_class(tibbleResult, "tbl_df")

    # Test with data.frame input
    dfInput <- as.data.frame(testSig)
    dfResult <- getConcordants(dfInput, "CP")
    expect_s3_class(dfResult, "tbl_df")

    # Test with S4Vectors::DataFrame input
    dataFrameInput <- S4Vectors::DataFrame(testSig)
    dataFrameResult <- getConcordants(dataFrameInput, "CP")
    expect_s4_class(dataFrameResult, "DFrame")
})

test_that("getConcordants numerical precision is maintained via VCR", {
    skip_on_cran()
    skip_on_bioc()
    skip_on_ci()

    vcr::local_cassette("getConcordants_precision_integration")

    testSig <- exampleSignature()[1L:3L, ]
    result <- getConcordants(testSig, ilincsLibrary = "CP")

    if (nrow(result) > 0L) {
        # Check similarity values are rounded to 8 decimal places
        similarityValues <- result[["similarity"]]
        for (val in similarityValues) {
            # Count decimal places
            decimalStr <- format(val, nsmall = 20L)
            decimalPart <- strsplit(decimalStr, "\\.")[[1L]][2L]
            if (!is.na(decimalPart)) {
                # Remove trailing zeros
                decimalPart <- gsub("0+$", "", decimalPart)
                expect_lte(nchar(decimalPart), 8L)
            }
        }

        # Check pValue precision (should be high precision)
        pValues <- result[["pValue"]]
        expect_true(all(is.numeric(pValues)))
        expect_true(all(pValues >= 0L & pValues <= 1L))
    }
})
