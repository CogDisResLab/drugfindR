# Comprehensive tests for consensusConcordants function and internal functions

# ==============================================================================
# TESTS FOR .validateConsensusConcordantsInput
# ==============================================================================

test_that(".validateConsensusConcordantsInput validates paired analysis requirements", {
    testConcordants <- getTestFixture("concordants_table",
        library = "CP",
        nRows = 100L, seed = .testSeed
    ) |>
        .processIlincsResponseSuccess("Up", "CP")

    # Paired analysis requires exactly 2 dataframes
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants), TRUE, 0.3, NULL),
        "Paired analysis requires two data frames"
    )

    # Paired analysis with 3 dataframes should error
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants, testConcordants, testConcordants), TRUE, 0.3, NULL),
        "Paired analysis requires two data frames"
    )

    # Valid paired analysis with 2 dataframes
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants, testConcordants), TRUE, 0.3, NULL)
    )
})

test_that(".validateConsensusConcordantsInput validates unpaired analysis requirements", {
    testConcordants <- getTestFixture("concordants_table",
        library = "CP",
        nRows = 100L, seed = .testSeed
    ) |>
        .processIlincsResponseSuccess("Up", "CP")

    # Unpaired analysis requires exactly 1 dataframe
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants, testConcordants), FALSE, 0.3, NULL),
        "Unpaired analysis requires only one dataframe"
    )

    # Unpaired analysis with no dataframes should error
    expect_error(
        .validateConsensusConcordantsInput(list(), FALSE, 0.3, NULL),
        "Unpaired analysis requires only one dataframe"
    )

    # Valid unpaired analysis with 1 dataframe
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0.3, NULL)
    )
})

test_that(".validateConsensusConcordantsInput validates cutoff parameter", {
    testConcordants <- getTestFixture("concordants_table",
        library = "CP",
        nRows = 100L, seed = .testSeed
    ) |>
        .processIlincsResponseSuccess("Up", "CP")

    # Non-numeric cutoff should error
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, "0.3", NULL),
        "cutoff must be a single numeric value"
    )

    # Multiple cutoff values should error
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, c(0.3, 0.5), NULL),
        "cutoff must be a single numeric value"
    )

    # Cutoff below 0 should error
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, -0.1, NULL),
        "cutoff must be between 0 and 1"
    )

    # Cutoff above 1 should error
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 1.5, NULL),
        "cutoff must be between 0 and 1"
    )

    # Valid cutoff values should work
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0.321, NULL)
    )
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0L, NULL)
    )
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 1L, NULL)
    )
})

test_that(".validateConsensusConcordantsInput validates cellLine parameter", {
    testConcordants <- getTestFixture("concordants_table",
        library = "CP",
        nRows = 100L, seed = .testSeed
    ) |>
        .processIlincsResponseSuccess("Up", "CP")

    # Non-character cellLine should error
    expect_error(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0.3, 123L),
        "cellLine must be a character vector or NULL"
    )

    # NULL cellLine should work
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0.3, NULL)
    )

    # Character vector cellLine should work
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0.3, "A375")
    )
    expect_silent(
        .validateConsensusConcordantsInput(list(testConcordants), FALSE, 0.3, c("A375", "PC3"))
    )
})

test_that(".validateConsensusConcordantsInput validates dataframe content", {
    # Empty dataframe should error
    emptyData <- getTestFixture("empty_ilincs_response") |>
        .processIlincsResponse(ilincsLibrary = "CP")

    expect_error(
        .validateConsensusConcordantsInput(list(emptyData), FALSE, 0.3, NULL),
        "All input dataframes must be non-empty"
    )

    # NULL dataframe should error
    expect_error(
        .validateConsensusConcordantsInput(list(NULL), FALSE, 0.3, NULL),
        "All input dataframes must be non-empty"
    )

    # Dataframe missing required columns should error
    incompleteData <- tibble::tibble(compound = "A", cellline = "A375")
    expect_error(
        .validateConsensusConcordantsInput(list(incompleteData), FALSE, 0.3, NULL),
        "Missing required columns"
    )
})

# ==============================================================================
# TESTS FOR .combineConcordantsData
# ==============================================================================

test_that(".combineConcordantsData combines single dataframe correctly", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .combineConcordantsData(list(testData))

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 10L)
    expect_identical(ncol(result), ncol(testData))
})

test_that(".combineConcordantsData combines multiple dataframes correctly", {
    testDataA <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Up", "CP")
    testDataB <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Down", "CP")

    result <- .combineConcordantsData(list(testDataA, testDataB))

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 20L) # 5 rows plus 3 rows
    expect_true(all(names(testDataA) %in% names(result)))
    expect_true(all(names(testDataB) %in% names(result)))
})

# ==============================================================================
# TESTS FOR .filterByCellLine
# ==============================================================================

test_that(".filterByCellLine returns original data when cellLine is NULL", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .filterByCellLine(testData, NULL)

    expect_identical(result, testData)
})

test_that(".filterByCellLine filters by single cell line correctly", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .filterByCellLine(testData, "HELA")

    expect_true(all(result[["cellline"]] == "HELA"))
    expect_lte(nrow(result), nrow(testData))
    expect_identical(nrow(result), 5L)
})

test_that(".filterByCellLine filters by multiple cell lines correctly", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .filterByCellLine(testData, c("HELA", "PC3"))

    expect_true(all(result[["cellline"]] %in% c("HELA", "PC3")))
    expect_false(any(result[["cellline"]] %in% c("A549", "A375")))
    expect_lte(nrow(result), nrow(testData))
    expect_identical(nrow(result), 6L)
})

test_that(".filterByCellLine returns empty result for non-existent cell lines", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")


    result <- .filterByCellLine(testData, "NONEXISTENT")

    expect_identical(nrow(result), 0L)
    expect_s3_class(result, "tbl_df")
})

# ==============================================================================
# TESTS FOR .applySimilarityCutoff
# ==============================================================================

test_that(".applySimilarityCutoff filters by absolute similarity correctly", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .applySimilarityCutoff(testData, 0.3)

    # Should keep |similarity| >= 0.3
    expect_identical(nrow(result), 6L) # 0.5, -0.8, 0.6, -0.3
    expect_true(all(abs(result[["similarity"]]) >= 0.3))
})

test_that(".applySimilarityCutoff handles zero cutoff", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .applySimilarityCutoff(testData, 0L)

    # Should keep all non-negative absolute values
    expect_identical(nrow(result), nrow(testData))
})

test_that(".applySimilarityCutoff handles high cutoff", {
    testData <- getTestFixture("concordants_table", seed = .testSeed) |>
        .processIlincsResponseSuccess("Any", "CP")

    result <- .applySimilarityCutoff(testData, 0.9)

    # Should return empty result
    expect_identical(nrow(result), 0L)
})

# ==============================================================================
# TESTS FOR .groupByTargetAndSelectMax
# ==============================================================================

test_that(".groupByTargetAndSelectMax selects maximum similarity per target", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, library = "KD") |>
        .processIlincsResponseSuccess("Any", "KD")

    result <- .groupByTargetAndSelectMax(testData)

    # Should have one entry per compound with max |similarity|
    expect_identical(nrow(result), 8L)

    # Check that we got the maximum similarities
    expect_identical(unique(sort(result[["treatment"]])), sort(result[["treatment"]]))
})

test_that(".groupByTargetAndSelectMax handles ties correctly", {
    testData <- tibble::tibble(
        signatureid = paste0("SIG", 1L:3L),
        treatment = c("A", "A", "A"),
        concentration = "10uM",
        time = "24h",
        cellline = c("A375", "PC3", "MCF7"),
        similarity = c(0.5, -0.5, 0.3), # Two with |similarity| = 0.5
        pValue = 0.01,
        sig_direction = "Up",
        sig_type = "Chemical Perturbagen"
    )

    result <- .groupByTargetAndSelectMax(testData)

    # Should keep both tied entries
    expect_identical(nrow(result), 2L)
    expect_true(all(abs(result[["similarity"]]) == 0.5))
})

# ==============================================================================
# TESTS FOR .selectAndOrderResults
# ==============================================================================

test_that(".selectAndOrderResults selects correct columns", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, library = "KD") |>
        .processIlincsResponseSuccess("Any", "KD")

    result <- .selectAndOrderResults(testData)

    expectedCols <- c(
        "signatureid", "treatment", "cellline", "time",
        "concentration", "sig_direction", "sig_type",
        "similarity", "pValue"
    )

    expect_named(result, expectedCols)
    expect_false("extra_column" %in% names(result))
})

test_that(".selectAndOrderResults orders by descending absolute similarity", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, library = "KD") |>
        .processIlincsResponseSuccess("Any", "KD")

    result <- .selectAndOrderResults(testData)

    # Should be ordered by descending |similarity|
    if (nrow(result) > 1L) {
        similarities <- abs(result[["similarity"]])
        expect_true(all(similarities[-length(similarities)] >= similarities[-1L]))
    }
})

# ==============================================================================
# TESTS FOR .applyTargetRenaming
# ==============================================================================

test_that(".applyTargetRenaming renames columns correctly for libraries", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, library = "KD") |>
        .processIlincsResponseSuccess("Any", "KD") |>
        .selectAndOrderResults()

    result <- .applyTargetRenaming(testData)

    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )

    expect_named(result, expectedCols)
})

# ==============================================================================
# TESTS FOR .processConsensusPipeline
# ==============================================================================

test_that(".processConsensusPipeline processes data through complete pipeline", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    result <- .processConsensusPipeline(testData, 0.3, "A375")

    # Should be filtered by cell line
    expect_true(all(result[["TargetCellLine"]] == "A375"))

    # Should be filtered by cutoff
    expect_true(all(abs(result[["Similarity"]]) >= 0.3))

    # Should have renamed columns
    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )

    expect_named(result, expectedCols)

    # Should be ordered by similarity
    if (nrow(result) > 1L) {
        similarities <- abs(result[["Similarity"]])
        expect_true(all(similarities[-length(similarities)] >= similarities[-1L]))
    }
})

test_that(".processConsensusPipeline handles NULL cellLine", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    result <- .processConsensusPipeline(testData, 0.3, NULL)

    expect_s3_class(result, "tbl_df")
    expect_true(all(abs(result[["Similarity"]]) >= 0.3))
})

# ==============================================================================
# TESTS FOR MAIN consensusConcordants FUNCTION
# ==============================================================================

test_that("consensusConcordants validates input correctly", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    # Test paired analysis validation
    expect_error(
        consensusConcordants(testData, paired = TRUE),
        "Paired analysis requires two data frames"
    )

    # Test unpaired analysis validation
    expect_error(
        consensusConcordants(testData, testData, paired = FALSE),
        "Unpaired analysis requires only one dataframe"
    )
})

test_that("consensusConcordants performs unpaired analysis correctly", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    result <- consensusConcordants(testData, cutoff = 0.321)

    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_true(all(abs(result[["Similarity"]]) >= 0.321))
    expect_named(result, expectedCols)
})

test_that("consensusConcordants performs paired analysis correctly", {
    testDataA <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Up", ilincsLibrary = "CP")
    testDataB <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Down", ilincsLibrary = "CP")

    result <- consensusConcordants(testDataA, testDataB, paired = TRUE)

    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
    expect_named(result, expectedCols)
})

test_that("consensusConcordants handles cell line filtering", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    result <- consensusConcordants(testData, cellLine = "A375")

    expect_true(all(result[["TargetCellLine"]] == "A375"))
    expect_identical(nrow(result), 20L)
})

test_that("consensusConcordants handles different cutoff values", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    resultLow <- consensusConcordants(testData, cutoff = 0.1)
    resultHigh <- consensusConcordants(testData, cutoff = 0.7)

    expect_lt(nrow(resultHigh), nrow(resultLow))
    expect_true(all(abs(resultHigh[["Similarity"]]) >= 0.7))
    expect_true(all(abs(resultLow[["Similarity"]]) >= 0.1))
    expect_identical(nrow(resultHigh), 27L)
    expect_identical(nrow(resultLow), 64L)
})

test_that("consensusConcordants handles empty results gracefully", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    result <- consensusConcordants(testData, cutoff = 1L)

    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )

    expect_named(result, expectedCols)
    expect_identical(nrow(result), 0L)
    expect_s3_class(result, "tbl_df")
})

test_that("consensusConcordants maintains data integrity", {
    testData <- getTestFixture("concordants_table", seed = .testSeed, nRows = 100L) |>
        .processIlincsResponseSuccess(sigDirection = "Any", ilincsLibrary = "CP")

    result <- consensusConcordants(testData)

    # Check data types
    expect_type(result[["TargetSignature"]], "character")
    expect_type(result[["Target"]], "character")
    expect_type(result[["Similarity"]], "double")
    expect_type(result[["pValue"]], "double")

    # Check no missing critical values
    expect_false(anyNA(result[["Target"]]))
    expect_false(any(is.infinite(result[["Similarity"]])))
})
