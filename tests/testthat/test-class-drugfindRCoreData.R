test_that("Constructor works with valid inputs", {
    sig <- exampleSignature()
    obj <- new("drugfindRCoreData",
        signature = sig,
        inputLibrary = "OE",
        outputLibrary = "CP",
        pairedAnalysis = FALSE,
        cellLines = c("MCF7", "PC3"),
        filterThresholdUp = 1.0,
        filterThresholdDown = -1.0,
        filteredSignature = sig,
        concordanceLimitUp = 0.5,
        concordanceLimitDown = -0.5,
        unfilteredConcordants = NULL,
        filteredConcordants = NULL
    )
    expect_s4_class(obj, "drugfindRCoreData")
    expect_equal(obj@signature, sig)
    expect_equal(obj@inputLibrary, "OE")
    expect_equal(obj@outputLibrary, "CP")
    expect_equal(obj@cellLines, c("MCF7", "PC3"))
    expect_false(obj@pairedAnalysis)
    expect_equal(obj@filterThresholdUp, 1.0)
    expect_equal(obj@filterThresholdDown, -1.0)
    expect_equal(obj@concordanceLimitUp, 0.5)
    expect_equal(obj@concordanceLimitDown, -0.5)
    expect_null(obj@unfilteredConcordants)
    expect_null(obj@filteredConcordants)
})

test_that("Default prototype values are set correctly", {
    obj <- new("drugfindRCoreData")
    expect_s4_class(obj, "drugfindRCoreData")
    expect_equal(obj@signature, tibble())
    expect_true(is.null(obj@inputLibrary))
    expect_true(is.null(obj@outputLibrary))
    expect_true(is.null(obj@cellLines))
    expect_equal(obj@filterThresholdUp, 0)
    expect_equal(obj@filterThresholdDown, 0)
    expect_equal(obj@concordanceLimitUp, 0.2)
    expect_equal(obj@concordanceLimitDown, -0.2)
    expect_true(obj@pairedAnalysis)
    expect_null(obj@filteredSignature)
    expect_null(obj@unfilteredConcordants)
    expect_null(obj@filteredConcordants)
})

test_that("Invalid library names fail validity check", {
    sig <- exampleSignature()
    obj <- new("drugfindRCoreData",
        signature = sig
    )

    obj@inputLibrary <- "XX"
    expect_match(
        validObject(obj, test = TRUE),
        "@inputLibrary must be one of 'OE', 'KD' or 'CP'"
    )
})

test_that("Concordance limits out of range fail validity", {
    sig <- exampleSignature()
    obj <- new("drugfindRCoreData",
        signature = sig
    )

    obj@concordanceLimitUp <- 2
    expect_match(
        validObject(obj, test = TRUE),
        "@concordanceLimitUp must be between -1 and 1"
    )

    obj@concordanceLimitUp <- 0.2
    obj@concordanceLimitDown <- -2
    expect_match(
        validObject(obj, test = TRUE),
        "@concordanceLimitDown must be between -1 and 1"
    )
})

test_that("Valid cell lines pass, invalid cell lines fail", {
    sig <- exampleSignature()
    valid <- new("drugfindRCoreData",
        signature = sig,
        inputLibrary = "OE",
        outputLibrary = "CP",
        pairedAnalysis = TRUE,
        cellLines = "MCF7",
        filterThresholdUp = 0.5,
        filterThresholdDown = -0.5,
        concordanceLimitUp = 0.5,
        concordanceLimitDown = -0.5
    )
    expect_true(validObject(valid, test = TRUE))

    invalid <- new("drugfindRCoreData",
        signature = sig
    )

    invalid@cellLines <- "MCCF"

    expect_match(validObject(invalid, test = TRUE), "@cellLines  must be one of the valid cell lines")
})
