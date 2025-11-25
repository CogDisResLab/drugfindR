# Comprehensive tests for filterSignature function and internal helpers

library(S4Vectors)

# ==============================================================================
# TESTS FOR INTERNAL VALIDATION FUNCTION
# ==============================================================================

test_that(
    ".validateFilterSignatureInput works correctly with valid inputs",
    {
        testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

        # Valid inputs should not error
        expect_silent(
            .validateFilterSignatureInput(testSignature, "any", 1.0, NULL)
        )
        expect_silent(
            .validateFilterSignatureInput(testSignature, "up", NULL, 0.1)
        )
        expect_silent(
            .validateFilterSignatureInput(
                testSignature, "down", c(1.0, 2.0),
                NULL
            )
        )

        # Test with different data frame types
        testDataFrame <- as.data.frame(testSignature)
        testS4DataFrame <- S4Vectors::DataFrame(testSignature)

        expect_silent(
            .validateFilterSignatureInput(testDataFrame, "any", 1.0, NULL)
        )
        expect_silent(
            .validateFilterSignatureInput(testS4DataFrame, "any", 1.0, NULL)
        )
    }
)

test_that(".validateFilterSignatureInput errors on invalid direction", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

    expect_error(
        .validateFilterSignatureInput(testSignature, "invalid", 1.0, NULL),
        "Direction must be one of 'up', 'down' or 'any'",
        fixed = TRUE
    )

    expect_error(
        .validateFilterSignatureInput(testSignature, "UP", 1.0, NULL),
        "Direction must be one of 'up', 'down' or 'any'",
        fixed = TRUE
    )

    expect_error(
        .validateFilterSignatureInput(testSignature, "", 1.0, NULL),
        "Direction must be one of 'up', 'down' or 'any'",
        fixed = TRUE
    )
})

test_that(".validateFilterSignatureInput errors on invalid threshold combinations", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

    # Both abs_threshold and pval_threshold specified should error
    expect_error(
        .validateFilterSignatureInput(testSignature, "any", 1.0, 0.1),
        "Only one of prop or threshold can be specified",
        fixed = TRUE
    )

    # Neither threshold specified should error
    expect_error(
        .validateFilterSignatureInput(testSignature, "any", NULL, NULL),
        "One of prop or threshold must be specified",
        fixed = TRUE
    )
})

test_that(".validateFilterSignatureInput errors on invalid abs_threshold", {
    testSignature <- getTestFixture("prepared_signature", seed = .testSeed)

    # Multiple thresholds should error (more than 2)
    expect_error(
        .validateFilterSignatureInput(testSignature, "any", c(1.0, 2.0, 3.0), NULL),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )

    # Empty threshold should error
    expect_error(
        .validateFilterSignatureInput(testSignature, "any", numeric(0L), NULL)
    )
})

test_that(".validateFilterSignatureInput errors on invalid threshold/prop combinations", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Both threshold and prop specified should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", 1.0, 0.1),
        "Only one of prop or threshold can be specified",
        fixed = TRUE
    )

    # Neither threshold nor prop specified should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", NULL, NULL),
        "One of prop or threshold must be specified",
        fixed = TRUE
    )
})

test_that(".validateFilterSignatureInput errors on invalid threshold length", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # More than two threshold values should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", c(1.0, 2.0, 3.0), NULL),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )

    # Empty threshold vector should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", numeric(0L), NULL),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )
})

test_that(".validateFilterSignatureInput errors on non-dataframe input", {
    # Non-dataframe input should error
    expect_error(
        .validateFilterSignatureInput("not_a_dataframe", "any", 1.0, NULL)
    )

    expect_error(
        .validateFilterSignatureInput(list(a = 1L, b = 2L), "any", 1.0, NULL)
    )

    expect_error(
        .validateFilterSignatureInput(matrix(1L:10L, nrow = 5L), "any", 1.0, NULL)
    )
})


test_that(".validateFilterSignatureInput errors on invalid threshold order", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Thresholds in wrong order should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", c(2.0, 1.0), NULL),
        "When two thresholds are specified, they must be in order \\(lower, higher\\)",
        fixed = FALSE
    )

    expect_error(
        .validateFilterSignatureInput(testSig, "any", c(0.5, -0.5), NULL),
        "When two thresholds are specified, they must be in order \\(lower, higher\\)",
        fixed = FALSE
    )
})

test_that(".validateFilterSignatureInput errors on invalid proportion values", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Proportion greater than 1 should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", NULL, 1.5),
        "Proportion must be between less than 0.5",
        fixed = TRUE
    )

    # Proportion equal to 0 should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", NULL, 0L),
        "Proportion must be between greater than 0",
        fixed = TRUE
    )

    # Proportion less than 0 should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", NULL, -0.1),
        "Proportion must be between greater than 0",
        fixed = TRUE
    )

    # Empty proportion vector should error
    expect_error(
        .validateFilterSignatureInput(testSig, "any", NULL, numeric(0L)),
        "Proportion must be specified as a single value",
        fixed = TRUE
    )
})

test_that(".validateFilterSignatureInput accepts valid proportion edge cases", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Proportion = 0.5 should be valid
    expect_silent(.validateFilterSignatureInput(testSig, "any", NULL, 0.5))
})

# ==============================================================================
# TESTS FOR THRESHOLD CALCULATION FUNCTIONS
# ==============================================================================

test_that("calculateSingleThreshold works correctly", {
    result <- .calculateSingleThreshold(1.5)

    expect_type(result, "list")
    expect_named(result, c("downThreshold", "upThreshold"))
    expect_identical(result[["downThreshold"]], -1.5)
    expect_identical(result[["upThreshold"]], 1.5)

    # Test with different values
    result2 <- .calculateSingleThreshold(0.8)
    expect_identical(result2[["downThreshold"]], -0.8)
    expect_identical(result2[["upThreshold"]], 0.8)

    # Test with zero
    result3 <- .calculateSingleThreshold(0.0)
    expect_identical(result3[["downThreshold"]], 0.0)
    expect_identical(result3[["upThreshold"]], 0.0)

    # Test with large value
    result4 <- .calculateSingleThreshold(100.0)
    expect_identical(result4[["downThreshold"]], -100.0)
    expect_identical(result4[["upThreshold"]], 100.0)
})

test_that(".calculateDoubleThreshold works correctly", {
    result <- .calculateDoubleThreshold(c(-2.0, 1.5))

    expect_type(result, "list")
    expect_named(result, c("downThreshold", "upThreshold"))
    expect_identical(result[["downThreshold"]], -2.0)
    expect_identical(result[["upThreshold"]], 1.5)

    # Test with different values
    result2 <- .calculateDoubleThreshold(c(-0.5, 2.5))
    expect_identical(result2[["downThreshold"]], -0.5)
    expect_identical(result2[["upThreshold"]], 2.5)

    # Test with identical values
    result3 <- .calculateDoubleThreshold(c(-1.0, 1.0))
    expect_identical(result3[["downThreshold"]], -1.0)
    expect_identical(result3[["upThreshold"]], 1.0)

    # Test that wrong order of thresholds errors
    expect_error(
        .calculateDoubleThreshold(c(1.0, -1.0)),
        "When two thresholds are specified, they must be in order (lower, higher)",
        fixed = TRUE
    )
})

test_that(".calculateAbsoluteThresholds dispatches correctly", {
    # Test single threshold dispatch
    result1 <- .calculateAbsoluteThresholds(1.0)
    expect_identical(result1[["downThreshold"]], -1.0)
    expect_identical(result1[["upThreshold"]], 1.0)

    # Test double threshold dispatch
    result2 <- .calculateAbsoluteThresholds(c(-1.5, 2.0))
    expect_identical(result2[["downThreshold"]], -1.5)
    expect_identical(result2[["upThreshold"]], 2.0)
})

test_that(".calculateAbsoluteThresholds errors on invalid input", {
    # Test invalid threshold (more than 2 values) should error
    expect_error(
        .calculateAbsoluteThresholds(c(1.0, 2.0, 3.0)),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )

    # Test invalid threshold (empty vector) should error
    expect_error(
        .calculateAbsoluteThresholds(numeric(0L)),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )

    # Test invalid threshold (4 values) should error
    expect_error(
        .calculateAbsoluteThresholds(1L:4L),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )
})

test_that(".calculateProportionalThreshold works correctly", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    calculateProp <- function(values, prop) {
        round(
            quantile(values, prop),
            2L
        )
    }

    # Test with 10% proportion
    result <- .calculateProportionalThreshold(testSig, 0.1)

    expect_type(result, "list")
    expect_named(result, c("downThreshold", "upThreshold"))

    # Check that thresholds are calculated correctly
    expectedDown <- calculateProp(testSig[["Value_LogDiffExp"]], 0.1)
    expectedUp <- calculateProp(testSig[["Value_LogDiffExp"]], 0.9)

    expect_identical(result[["downThreshold"]], expectedDown)
    expect_identical(result[["upThreshold"]], expectedUp)

    # Test with different proportion
    result2 <- .calculateProportionalThreshold(testSig, 0.2)
    expectedDown2 <- calculateProp(testSig[["Value_LogDiffExp"]], 0.2)
    expectedUp2 <- calculateProp(testSig[["Value_LogDiffExp"]], 0.8)

    expect_identical(result2[["downThreshold"]], expectedDown2)
    expect_identical(result2[["upThreshold"]], expectedUp2)

    # Test with extreme proportion
    result3 <- .calculateProportionalThreshold(testSig, 0.05)
    expectedDown3 <- calculateProp(testSig[["Value_LogDiffExp"]], 0.05)
    expectedUp3 <- calculateProp(testSig[["Value_LogDiffExp"]], 0.95)

    expect_identical(result3[["downThreshold"]], expectedDown3)
    expect_identical(result3[["upThreshold"]], expectedUp3)

    # Test with 50% proportion
    result4 <- .calculateProportionalThreshold(testSig, 0.5)
    expectedMedian <- calculateProp(testSig[["Value_LogDiffExp"]], 0.5)

    expect_identical(result4[["downThreshold"]], expectedMedian)
    expect_identical(result4[["upThreshold"]], expectedMedian)
})

# ==============================================================================
# TESTS FOR DIRECTION FILTERING FUNCTION
# ==============================================================================

test_that(".applyDirectionFilter works correctly for 'up' direction", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)
    thresholds <- list(downThreshold = -1.5, upThreshold = 1.5)

    result <- .applyDirectionFilter(testSig, thresholds, "up")

    expect_true(all(result[["Value_LogDiffExp"]] >= 1.5))
    expect_identical(nrow(result), 21L) # Values: 2.0, 3.0, 4.0

    # Check specific values
    expect_setequal(
        result[["Value_LogDiffExp"]],
        .applyDirectionFilter(
            testSig,
            thresholds,
            "up"
        )[["Value_LogDiffExp"]]
    )
    expect_setequal(
        result[["Name_GeneSymbol"]],
        .applyDirectionFilter(
            testSig,
            thresholds,
            "up"
        )[["Name_GeneSymbol"]]
    )
})

test_that(".applyDirectionFilter works correctly for 'down' direction", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)
    thresholds <- list(downThreshold = -1.5, upThreshold = 1.5)

    result <- .applyDirectionFilter(testSig, thresholds, "down")

    expect_true(all(result[["Value_LogDiffExp"]] <= -1.5))
    expect_identical(nrow(result), 21L) # Values: -3.0, -2.0

    # Check specific values
    expect_setequal(
        result[["Value_LogDiffExp"]],
        .applyDirectionFilter(
            testSig,
            thresholds,
            "down"
        )[["Value_LogDiffExp"]]
    )
    expect_setequal(
        result[["Name_GeneSymbol"]],
        .applyDirectionFilter(
            testSig,
            thresholds,
            "down"
        )[["Name_GeneSymbol"]]
    )
})

test_that(".applyDirectionFilter works correctly for 'any' direction", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)
    thresholds <- list(downThreshold = -1.5, upThreshold = 1.5)

    result <- .applyDirectionFilter(testSig, thresholds, "any")

    expect_true(all(
        result[["Value_LogDiffExp"]] >= 1.5 |
            result[["Value_LogDiffExp"]] <= -1.5
    ))
    expect_identical(nrow(result), 42L) # Values: -3.0, -2.0, 2.0, 3.0, 4.0

    # Check specific values
    expect_setequal(
        result[["Value_LogDiffExp"]],
        .applyDirectionFilter(testSig, thresholds, "any")[["Value_LogDiffExp"]]
    )
    expect_setequal(
        result[["Name_GeneSymbol"]],
        .applyDirectionFilter(testSig, thresholds, "any")[["Name_GeneSymbol"]]
    )
})

test_that(".applyDirectionFilter works with edge case thresholds", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)

    # Test with very high thresholds (should return empty)
    highThresholds <- list(downThreshold = -10.0, upThreshold = 10.0)
    resultEmpty <- .applyDirectionFilter(testSig, highThresholds, "any")
    expect_identical(nrow(resultEmpty), 0L)

    # Test with zero thresholds
    zeroThresholds <- list(downThreshold = 0.0, upThreshold = 0.0)
    resultZeroAny <- .applyDirectionFilter(testSig, zeroThresholds, "any")
    resultZeroUp <- .applyDirectionFilter(testSig, zeroThresholds, "up")
    resultZeroDown <- .applyDirectionFilter(testSig, zeroThresholds, "down")

    expect_identical(nrow(resultZeroAny), 100L) # All values
    expect_identical(nrow(resultZeroUp), 55L) # Values >= 0: 0.0, 1.0, 2.0, 3.0, 4.0
    expect_identical(nrow(resultZeroDown), 45L) # Values <= 0: -3.0, -2.0, -1.0, -0.5, 0.0
})

test_that(".applyDirectionFilter maintains data structure", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)
    thresholds <- list(downThreshold = -1.5, upThreshold = 1.5)

    result <- .applyDirectionFilter(testSig, thresholds, "any")

    # Should maintain all columns
    expect_identical(colnames(result), colnames(testSig))

    # Should maintain column types
    expect_identical(
        vapply(result, class, character(1L)),
        vapply(testSig, class, character(1L))
    )
})

# ==============================================================================
# TESTS FOR MAIN FUNCTION INPUT VALIDATION
# ==============================================================================

test_that("filterSignature errors on invalid input combinations", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Not specifying threshold and prop causes error
    expect_error(
        filterSignature(testSig),
        "One of prop or threshold must be specified",
        fixed = TRUE
    )

    # Specifying both threshold and prop causes error
    expect_error(
        filterSignature(testSig, threshold = 0.1, prop = 0.1),
        "Only one of prop or threshold can be specified",
        fixed = TRUE
    )

    # Invalid signature direction causes error
    expect_error(
        filterSignature(testSig, direction = "invalid", threshold = 1.0),
        "Direction must be one of 'up', 'down' or 'any'",
        fixed = TRUE
    )

    # More than two threshold values causes error
    expect_error(
        filterSignature(testSig, threshold = c(0.0, 0.1, 0.2)),
        "Threshold must be specified as one or two values",
        fixed = TRUE
    )
})

test_that("filterSignature handles empty signatures", {
    emptySig <- createEmptySignature()

    result <- filterSignature(emptySig, threshold = 0.0)
    expect_identical(nrow(result), 0L)
})

# ==============================================================================
# TESTS FOR THRESHOLD-BASED FILTERING
# ==============================================================================

test_that("filterSignature works with single threshold value", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)

    # Test with zero threshold
    resultZero <- filterSignature(testSig, threshold = 0.0)
    expect_identical(nrow(resultZero), 100L) # All values should pass

    # Test with non-zero threshold
    resultNonzero <- filterSignature(testSig, threshold = 1.5)
    expect_true(all(abs(resultNonzero[["Value_LogDiffExp"]]) >= 1.5))
    expect_identical(nrow(resultNonzero), 42L) # Values: -3, -2, -1, 1, 2, 3, 4

    # Test with high threshold
    resultHigh <- filterSignature(testSig, threshold = 2.5)
    expect_true(all(abs(resultHigh[["Value_LogDiffExp"]]) >= 2.5))
    expect_identical(nrow(resultHigh), 12L) # Values: -3, 3, 4
})

test_that("filterSignature works with single threshold and direction", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test up-regulated genes only
    resultUp <- filterSignature(testSig, threshold = 1.5, direction = "up")
    expect_true(all(resultUp[["Value_LogDiffExp"]] >= 1.5))
    expect_identical(nrow(resultUp), 1L) # Values: 1, 2, 3, 4

    # Test down-regulated genes only
    resultDown <- filterSignature(testSig, threshold = 1.5, direction = "down")
    expect_true(all(resultDown[["Value_LogDiffExp"]] <= -1.5))
    expect_identical(nrow(resultDown), 2L) # Values: -3, -2, -1

    # Test both directions explicitly
    resultAny <- filterSignature(testSig, threshold = 1.5, direction = "any")
    expect_true(all(abs(resultAny[["Value_LogDiffExp"]]) >= 1.5))
    expect_identical(nrow(resultAny), 3L) # Values: -3, -2, -1, 1, 2, 3, 4
})

test_that("filterSignature works with double threshold values", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test asymmetric thresholds
    result <- filterSignature(testSig, threshold = c(-0.75, 1.5))

    # Check up-regulated genes
    upGenes <- result[result[["Value_LogDiffExp"]] >= 0L, ]
    expect_true(all(upGenes[["Value_LogDiffExp"]] >= 1.5))

    # Check down-regulated genes
    downGenes <- result[result[["Value_LogDiffExp"]] <= 0L, ]
    expect_true(all(downGenes[["Value_LogDiffExp"]] <= -0.75))

    expect_identical(nrow(result), 3L) # Values: -3, -2, -1, 2, 3, 4
})

test_that("filterSignature works with double threshold and direction", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test up-regulated genes with asymmetric thresholds
    resultUp <- filterSignature(testSig, threshold = c(-2.5, 1.5), direction = "up")
    expect_true(all(resultUp[["Value_LogDiffExp"]] >= 1.5))
    expect_identical(nrow(resultUp), 1L) # Values: 2, 3, 4

    # Test down-regulated genes with asymmetric thresholds
    resultDown <- filterSignature(testSig, threshold = c(-1.5, 2.5), direction = "down")
    expect_true(all(resultDown[["Value_LogDiffExp"]] <= -1.5))
    expect_identical(nrow(resultDown), 2L) # Values: -3, -2
})

# ==============================================================================
# TESTS FOR PROPORTION-BASED FILTERING
# ==============================================================================

test_that("filterSignature works with proportion filtering", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)

    # Test with 100% proportion (should include all)
    resultAll <- filterSignature(testSig, prop = 0.5)
    expect_identical(nrow(resultAll), 100L)

    # Test with specific proportion
    result20 <- filterSignature(testSig, prop = 0.2)
    expectedDown <- quantile(testSig[["Value_LogDiffExp"]], 0.2)
    expectedUp <- quantile(testSig[["Value_LogDiffExp"]], 0.8)

    upGenes <- result20[result20[["Value_LogDiffExp"]] >= 0L, ]
    downGenes <- result20[result20[["Value_LogDiffExp"]] <= 0L, ]

    expect_true(all(upGenes[["Value_LogDiffExp"]] >= expectedUp))
    expect_true(all(downGenes[["Value_LogDiffExp"]] <= expectedDown))
})

test_that("filterSignature works with proportion and direction", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed, nGenes = 100L)

    calculateProp <- function(values, prop) {
        round(
            quantile(values, prop),
            2L
        )
    }

    # Test up-regulated genes only with proportion
    resultUp <- filterSignature(testSig, prop = 0.1, direction = "up")
    expectedUp <- calculateProp(testSig[["Value_LogDiffExp"]], 0.9)

    expect_true(all(resultUp[["Value_LogDiffExp"]] >= expectedUp))
    expect_identical(nrow(resultUp), 10L) # Top 10% = 1 gene

    # Test down-regulated genes only with proportion
    resultDown <- filterSignature(testSig, prop = 0.1, direction = "down")
    expectedDown <- calculateProp(testSig[["Value_LogDiffExp"]], 0.1)

    expect_true(all(resultDown[["Value_LogDiffExp"]] <= expectedDown))
    expect_identical(nrow(resultDown), 10L) # Bottom 10% = 1 gene
})

# ==============================================================================
# TESTS FOR DATA FRAME COMPATIBILITY
# ==============================================================================

test_that("filterSignature maintains compatibility with different data frame types", {
    # Create test data in different formats
    testTypeTibble <- getTestFixture("prepared_signature", seed = .testSeed)
    testTypeDf <- as.data.frame(testTypeTibble)
    testTypeDataFrame <- S4Vectors::DataFrame(testTypeTibble)

    # Apply same filtering to all
    resultTibble <- filterSignature(testTypeTibble, threshold = 1.0)
    resultDf <- filterSignature(testTypeDf, threshold = 1.0)
    resultDataFrame <- filterSignature(testTypeDataFrame, threshold = 1.0)

    # All should return tibbles
    expect_s3_class(resultTibble, "tbl_df")
    expect_s3_class(resultDf, "data.frame")
    expect_s4_class(resultDataFrame, "DFrame")

    # Results should be identical in content
    expect_identical(as.data.frame(resultTibble), as.data.frame(resultDf))
    expect_identical(as.data.frame(resultTibble), as.data.frame(resultDataFrame))
    expect_identical(as.data.frame(resultDf), as.data.frame(resultDataFrame))

    # Check that all have same genes
    expect_setequal(resultTibble[["Name_GeneSymbol"]], resultDf[["Name_GeneSymbol"]])
    expect_setequal(resultTibble[["Name_GeneSymbol"]], resultDataFrame[["Name_GeneSymbol"]])
})

test_that("filterSignature preserves data frame structure", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    result <- filterSignature(testSig, threshold = 1.0)

    # Should maintain column structure
    expect_identical(colnames(result), colnames(testSig))
    expect_identical(
        vapply(result, class, character(1L)),
        vapply(testSig, class, character(1L))
    )

    # Should maintain row structure integrity
    expect_true(all(c(
        "Name_GeneSymbol",
        "Value_LogDiffExp",
        "Significance_pvalue"
    ) %in% colnames(result)))
})

# ==============================================================================
# TESTS FOR EDGE CASES AND ERROR HANDLING
# ==============================================================================

test_that("filterSignature handles edge cases correctly", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test with very high threshold (should return empty result)
    resultEmpty <- filterSignature(testSig, threshold = 10.0)
    expect_identical(nrow(resultEmpty), 0L)
    expect_s3_class(resultEmpty, "tbl_df")
    expect_identical(colnames(resultEmpty), colnames(testSig))

    # Test with zero threshold
    resultZero <- filterSignature(testSig, threshold = 0.0)
    expect_identical(nrow(resultZero), nrow(testSig)) # All rows should pass

    # Test with very small proportion
    resultSmallProp <- filterSignature(testSig, prop = 0.01)
    expect_gte(nrow(resultSmallProp), 0L)
    expect_s3_class(resultSmallProp, "tbl_df")

    # Test with maximum proportion
    resultMaxProp <- filterSignature(testSig, prop = 0.5)
    expect_identical(nrow(resultMaxProp), nrow(testSig)) # Should return all data
})

# ==============================================================================
# COMPREHENSIVE INTEGRATION TESTS
# ==============================================================================

test_that("filterSignature refactored version maintains backward compatibility", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test all combinations of parameters
    testCases <- list(
        list(threshold = 1.5, direction = "any"),
        list(threshold = 1.5, direction = "up"),
        list(threshold = 1.5, direction = "down"),
        list(threshold = c(-2.0, 1.0), direction = "any"),
        list(threshold = c(-2.0, 1.0), direction = "up"),
        list(threshold = c(-2.0, 1.0), direction = "down"),
        list(prop = 0.1, direction = "any"),
        list(prop = 0.1, direction = "up"),
        list(prop = 0.1, direction = "down"),
        list(prop = 0.2, direction = "any")
    )

    for (case in testCases) {
        result <- do.call(filterSignature, c(list(signature = testSig), case))

        # Should always return a tibble
        expect_s3_class(result, "tbl_df")

        # Should maintain column structure
        expect_identical(colnames(result), colnames(testSig))

        # Should have appropriate number of rows
        expect_gte(nrow(result), 0L)
        expect_lte(nrow(result), nrow(testSig))

        # Should maintain data integrity
        if (nrow(result) > 0L) {
            expect_true(all(result[["signatureID"]] == "LINCSID_00000"))
            expect_true(all(result[["Name_GeneSymbol"]] %in% testSig[["Name_GeneSymbol"]]))
        }
    }
})
