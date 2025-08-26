# r Comprehensive tests for getSignature function and internal helpers

library(tibble)
library(httr2)
library(vcr)
library(dplyr)

# ==============================================================================
# TESTS FOR INTERNAL VALIDATION FUNCTION
# ==============================================================================

test_that(".validateGetSignatureInput works correctly with valid inputs", {
    # Valid signature IDs should not error
    expect_silent(.validateGetSignatureInput(lincsKdId()))
    expect_silent(.validateGetSignatureInput(lincsOeId()))
    expect_silent(.validateGetSignatureInput(lincsCpId()))
})

test_that(".validateGetSignatureInput errors on invalid sigId type", {
    # Non-character input should error
    expect_error(
        .validateGetSignatureInput(123L), # nolint: implicit_integer_linter.
        "sigId must be a character string"
    )

    expect_error(
        .validateGetSignatureInput(TRUE),
        "sigId must be a character string"
    )

    expect_error(
        .validateGetSignatureInput(list(lincsKdId())),
        "sigId must be a character string"
    )
})

test_that(".validateGetSignatureInput errors on invalid sigId length", {
    # Multiple values should error
    expect_error(
        .validateGetSignatureInput(c(lincsKdId(), "LINCSOE_123")),
        "sigId must be a single character string"
    )

    expect_error(
        .validateGetSignatureInput(character(0L)),
        "sigId must be a single character string"
    )
})

test_that(".validateGetSignatureInput errors on empty or whitespace sigId", {
    # Empty string should error
    expect_error(
        .validateGetSignatureInput(""),
        "sigId cannot be empty or consist only of whitespace"
    )

    # Whitespace-only string should error
    expect_error(
        .validateGetSignatureInput("   "),
        "sigId cannot be empty or consist only of whitespace"
    )

    expect_error(
        .validateGetSignatureInput("\t\n  "),
        "sigId cannot be empty or consist only of whitespace"
    )
})

test_that(".validateGetSignatureInput errors on nonexistent signature", {
    # Signature not in metadata should error
    expect_error(
        .validateGetSignatureInput("NONEXISTENT_SIG"),
        "Signature ID 'NONEXISTENT_SIG' not found in metadata tables"
    )

    expect_error(
        .validateGetSignatureInput("INVALID_FORMAT"),
        "Signature ID 'INVALID_FORMAT' not found in metadata tables"
    )
})

# ==============================================================================
# TESTS FOR SIGNATURE EXISTENCE VALIDATION
# ==============================================================================

test_that(".isValidSignatureId correctly identifies valid signatures", {
    # Test with known valid signatures from each library
    expect_true(.isValidSignatureId(lincsKdId()))

    # Test with first signatures from each metadata table
    if (nrow(cpMetadata) > 0L) {
        firstCpSignature <- cpMetadata[["SourceSignature"]][10L]
        expect_true(.isValidSignatureId(firstCpSignature))
    }

    if (nrow(kdMetadata) > 0L) {
        firstKdSignature <- kdMetadata[["SourceSignature"]][10L]
        expect_true(.isValidSignatureId(firstKdSignature))
    }

    if (nrow(oeMetadata) > 0L) {
        firstOeSignature <- oeMetadata[["SourceSignature"]][10L]
        expect_true(.isValidSignatureId(firstOeSignature))
    }
})

test_that(".isValidSignatureId correctly identifies invalid signatures", {
    # Test with clearly invalid signatures
    expect_false(.isValidSignatureId("NONEXISTENT_SIG"))
    expect_false(.isValidSignatureId("INVALID_FORMAT"))
    expect_false(.isValidSignatureId(""))
    expect_false(.isValidSignatureId("NOT_A_REAL_SIGNATURE"))
})

# ==============================================================================
# TESTS FOR L1000 DETECTION FUNCTION
# ==============================================================================

test_that(".detectL1000Status correctly identifies L1000 signatures", {
    # Test with known L1000 signatures
    expect_true(.detectL1000Status(lincsKdId()))

    # Test pattern matching for L1000 signatures
    testL1000Patterns <- c(lincsKdId(), lincsCpId(), lincsOeId())
    for (pattern in testL1000Patterns) {
        if (.isValidSignatureId(pattern)) {
            expect_true(.detectL1000Status(pattern))
        }
    }
})

test_that(".detectL1000Status correctly identifies non-L1000 signatures", {
    # Test with signatures that don't match L1000 pattern
    expect_false(.detectL1000Status("GEO_UP_1"))
    expect_false(.detectL1000Status("CUSTOM_SIG_123"))
    expect_false(.detectL1000Status("NOT_L1000_FORMAT"))

    # Test with invalid signatures (should return FALSE)
    expect_false(.detectL1000Status("NONEXISTENT_SIG"))
    expect_false(.detectL1000Status("LINCSKD_NOTANUMBER"))
})

test_that(".detectL1000Status validates signature existence", {
    # Should return FALSE for patterns that look like L1000 but don't exist
    expect_false(.detectL1000Status("LINCSKD_999999"))
    expect_false(.detectL1000Status("LINCSOE_999999"))
    expect_false(.detectL1000Status("LINCSCP_999999"))
})

# ==============================================================================
# TESTS FOR HTTP REQUEST CREATION
# ==============================================================================

test_that(".createSignatureRequest creates valid httr2 request object", {
    # Test with valid signature ID
    request <- .createSignatureRequest(lincsKdId())

    # Should return an httr2_request object
    expect_s3_class(request, "httr2_request")

    # Check URL construction
    requestUrl <- req_get_url(request)
    parsedRequestUrl <- url_parse(requestUrl)
    expect_identical(parsedRequestUrl[["scheme"]], "https")
    expect_identical(parsedRequestUrl[["hostname"]], "www.ilincs.org")
    expect_identical(parsedRequestUrl[["path"]], "/api/ilincsR/downloadSignature")

    # Check query parameters
    queryParameters <- parsedRequestUrl[["query"]]
    expect_identical(queryParameters[["sigID"]], lincsKdId())
    expect_identical(queryParameters[["noOfTopGenes"]], "Inf")

    # Check method
    expect_identical(req_get_method(request), "GET")

    # Check user agent
    userAgent <- request[["options"]][["useragent"]]
    expect_true(grepl("drugfindR", userAgent))
    expect_true(grepl("github.com/CogDisResLab/drugfindR", userAgent))
})

test_that(".createSignatureRequest works with different signature IDs", {
    # Test with different types of valid signature IDs
    testSignatureIds <- c(lincsKdId(), lincsOeId(), lincsCpId())

    for (signatureId in testSignatureIds) {
        if (.isValidSignatureId(signatureId)) {
            request <- .createSignatureRequest(signatureId)
            requestUrl <- req_get_url(request)
            parsedRequestUrl <- url_parse(requestUrl)
            queryParameters <- parsedRequestUrl[["query"]]

            expect_s3_class(request, "httr2_request")
            expect_identical(queryParameters[["sigID"]], signatureId)
        }
    }
})

# ==============================================================================
# TESTS FOR HTTP REQUEST EXECUTION
# ==============================================================================

test_that(".executeSignatureRequest configures request properly", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        request <- .createSignatureRequest(lincsKdId())
        response <- .executeSignatureRequest(request)
        expect_s3_class(response, "httr2_response")
    })
})

test_that(".executeSignatureRequest handles verbose option", {
    # TODO: Ensure this test works
    skip()
    successResponse <- mockSignatureResponse("LINCSKD_28", 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        request <- .createSignatureRequest("LINCSKD_28")

        # Test verbose = FALSE (default)
        response1 <- .executeSignatureRequest(request)
        response1Json <- resp_body_json(response1)
        expect_null(response1Json[["request"]][["options"]][["verbose"]])
        expect_s3_class(response1, "httr2_response")

        # Test verbose = FALSE
        response2 <- .executeSignatureRequest(request, verbose = FALSE)
        response2Json <- resp_body_json(response2)
        expect_null(response2Json[["request"]][["options"]][["verbose"]])
        expect_s3_class(response1, "httr2_response")

        # Test verbose = TRUE
        response3 <- .executeSignatureRequest(request, verbose = TRUE)
        response3Json <- resp_body_json(response3)
        expect_true(response3Json[["request"]][["options"]][["verbose"]])
        expect_s3_class(response3, "httr2_response")
    })
})

test_that(".executeSignatureRequest returns error responses without raising", {
    errorResponse <- mockSignatureResponse("", 400L, seed = .testSeed)

    with_mocked_responses(errorResponse, {
        request <- .createSignatureRequest("")
        response <- .executeSignatureRequest(request)
        expect_s3_class(response, "httr2_response")
        expect_identical(httr2::resp_status(response), 400L)
    })
})

# ==============================================================================
# TESTS FOR SUCCESSFUL RESPONSE PROCESSING
# ==============================================================================

test_that(".processSuccessfulResponse handles valid response correctly", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        request <- .createSignatureRequest(lincsKdId())
        response <- .executeSignatureRequest(request)
        result <- .processSuccessfulResponse(response)

        # Should return a tibble
        expect_s3_class(result, "tbl_df")

        # Check column names and structure
        expectedColumns <- signatureColumns()
        expect_named(result, expectedColumns)

        # Check data types
        expect_type(result[["signatureID"]], "character")
        expect_type(result[["ID_geneid"]], "character")
        expect_type(result[["Name_GeneSymbol"]], "character")
        expect_type(result[["Value_LogDiffExp"]], "double")
        expect_type(result[["Significance_pvalue"]], "double")
        expect_type(result[["is_L1000"]], "logical")
    })
})

test_that(".processSuccessfulResponse adds correct metadata", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        request <- .createSignatureRequest(lincsKdId())
        response <- .executeSignatureRequest(request)
        result <- .processSuccessfulResponse(response)

        # Check signatureID is correctly added
        expect_true(all(result[["signatureID"]] == lincsKdId()))

        # Check L1000 status is correctly detected
        expect_true(all(result[["is_L1000"]] == TRUE)) # Signature is L1000
    })
})

test_that(".processSuccessfulResponse rounds numerical values correctly", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        request <- .createSignatureRequest(lincsKdId())
        response <- .executeSignatureRequest(request)
        result <- .processSuccessfulResponse(response)

        # Check that values are rounded to 12 decimal places
        logFcDigits <- vapply(result[["Value_LogDiffExp"]], function(x) {
            decimalPart <- x - floor(x)
            if (decimalPart == 0L) {
                return(0L)
            }
            nchar(gsub("0+$", "", format(decimalPart, scientific = FALSE))) - 2L
        }, FUN.VALUE = numeric(1L))
        expect_true(all(logFcDigits <= 12L))

        pValueDigits <- vapply(result[["Significance_pvalue"]], function(x) {
            decimalPart <- x - floor(x)
            if (decimalPart == 0L) {
                return(0L)
            }
            nchar(gsub("0+$", "", format(decimalPart, scientific = FALSE))) - 2L
        }, FUN.VALUE = numeric(1L))
        expect_true(all(pValueDigits <= 12L))
    })
})

# ==============================================================================
# TESTS FOR ERROR RESPONSE PROCESSING
# ==============================================================================

test_that(".processSignatureResponseError handles 400 errors correctly", {
    errorResponse <- mockSignatureResponse("", 400L, seed = .testSeed)

    with_mocked_responses(errorResponse, {
        request <- .createSignatureRequest("LINCSKD_28")
        response <- .executeSignatureRequest(request)

        expect_error(
            .processSignatureResponseError(response),
            "iLINCS API request failed \\(Status 400\\): Bad Request"
        )
    })
})

test_that(".processSignatureResponseError handles 500 errors correctly", {
    errorResponse <- mockSignatureResponse("", 500L, seed = .testSeed)

    with_mocked_responses(errorResponse, {
        request <- .createSignatureRequest("LINCSKD_28")
        response <- .executeSignatureRequest(request)

        expect_error(
            .processSignatureResponseError(response),
            "iLINCS API request failed \\(Status 500\\): Internal Server Error"
        )
    })
})

# ==============================================================================
# TESTS FOR RESPONSE PROCESSING DISPATCHER
# ==============================================================================

test_that(".processSignatureResponse dispatcher works correctly for success case", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        request <- .createSignatureRequest(lincsKdId())
        response <- .executeSignatureRequest(request)
        result <- .processSignatureResponse(response)

        # Should return a tibble with expected structure
        expect_s3_class(result, "tbl_df")
        expect_named(result, signatureColumns())
    })
})

test_that(".processSignatureResponse dispatcher works correctly for error case", {
    errorResponse <- mockSignatureResponse("", 404L, seed = .testSeed)

    with_mocked_responses(errorResponse, {
        request <- .createSignatureRequest(lincsKdId())
        response <- .executeSignatureRequest(request)

        expect_error(
            .processSignatureResponse(response),
            "iLINCS API request failed"
        )
    })
})

# ==============================================================================
# VCR INTEGRATION TESTS - API-LEVEL CHECKS
# ==============================================================================

test_that("getSignature works end-to-end with valid KD signature via VCR", {
    vcr::local_cassette("getSignature_valid_KD")
    result <- getSignature(lincsKdId())

    # Check basic structure
    expect_s3_class(result, "tbl_df")
    expect_named(result, signatureColumns())

    # Check data completeness
    expect_identical(nrow(result), 978L)
    expect_false(any(purrr::flatten_lgl(purrr::map(result, is.na))))

    # Check L1000 status
    expect_true(all(result[["is_L1000"]] == TRUE))

    # Check signature ID consistency
    expect_true(all(result[["signatureID"]] == lincsKdId()))

    # Check data types
    expect_type(result[["ID_geneid"]], "character")
    expect_type(result[["Name_GeneSymbol"]], "character")
    expect_type(result[["Value_LogDiffExp"]], "double")
    expect_type(result[["Significance_pvalue"]], "double")
    expect_type(result[["is_L1000"]], "logical")
})

# ==============================================================================
# TESTS FOR EDGE CASES AND ERROR HANDLING
# ==============================================================================

test_that("getSignature maintains numerical precision", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        result <- getSignature(lincsKdId())

        # Values should be rounded to 12 decimal places max
        logFcValues <- result[["Value_LogDiffExp"]]
        pValueValues <- result[["Significance_pvalue"]]

        # Check that rounding is applied consistently
        expect_true(all(logFcValues == round(logFcValues, 12L)))
        expect_true(all(pValueValues == round(pValueValues, 12L)))
    })
})

test_that("getSignature handles different input validation scenarios", {
    # Test various invalid inputs
    expect_error(getSignature(123L), "sigId must be a character string")
    expect_error(getSignature(c("A", "B")), "sigId must be a single character string")
    expect_error(getSignature(""), "sigId cannot be empty")
    expect_error(getSignature("   "), "sigId cannot be empty")
})

test_that("getSignature output structure is consistent", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        result <- getSignature(lincsKdId())

        # Check that output always has the same column structure
        expectedColumns <- signatureColumns()
        expect_identical(colnames(result), expectedColumns)

        # Check that columns are in the expected order
        expect_identical(colnames(result)[1L], "signatureID")
        expect_identical(colnames(result)[6L], "is_L1000")
    })
})

# ==============================================================================
# PERFORMANCE AND INTEGRATION TESTS
# ==============================================================================

test_that("getSignature works with realistic data volumes", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        result <- getSignature(lincsKdId())

        # Should handle standard L1000 gene set size
        expect_gte(nrow(result), 1L)
        expect_lte(nrow(result), 978L) # Reasonable upper bound

        # All rows should have complete data
        expect_false(any(is.na(result[["signatureID"]])))
        expect_false(any(is.na(result[["ID_geneid"]])))
        expect_false(any(is.na(result[["Name_GeneSymbol"]])))
    })
})

test_that("getSignature is compatible with downstream functions", {
    successResponse <- mockSignatureResponse(lincsKdId(), 200L, seed = .testSeed)

    with_mocked_responses(successResponse, {
        signature <- getSignature(lincsKdId())

        # Should be compatible with filterSignature (after removing is_L1000 column)
        filterableSignature <- dplyr::select(signature, -"is_L1000")
        expect_silent(stopIfInvalidSignature(filterableSignature))

        # Should have the correct structure for downstream analysis
        expect_true("signatureID" %in% colnames(signature))
        expect_true("ID_geneid" %in% colnames(signature))
        expect_true("Name_GeneSymbol" %in% colnames(signature))
        expect_true("Value_LogDiffExp" %in% colnames(signature))
        expect_true("Significance_pvalue" %in% colnames(signature))
        expect_true("is_L1000" %in% colnames(signature))
    })
})
