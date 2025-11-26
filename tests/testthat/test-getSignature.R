# r Comprehensive tests for getSignature function and internal helpers

library(httr2)

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
        .validateGetSignatureInput(c(lincsKdId(), lincsCpId())),
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
    # nolint start: nonportable_path_linter, absolute_path_linter
    expect_identical(parsedRequestUrl[["path"]], "/api/ilincsR/downloadSignature")
    # nolint end

    # Check query parameters
    queryParameters <- parsedRequestUrl[["query"]]
    expect_identical(queryParameters[["sigID"]], lincsKdId())
    expect_identical(queryParameters[["noOfTopGenes"]], "Inf")

    # Check method
    expect_identical(req_get_method(request), "GET")

    # Check user agent
    userAgent <- request[["options"]][["useragent"]]
    expect_true(grepl("github.com/CogDisResLab/drugfindR", userAgent)) # nolint: nonportable_path_linter.
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
