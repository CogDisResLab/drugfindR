# Comprehensive tests for consensusConcordants function and internal functions

library(tibble)
library(dplyr)

# ==============================================================================
# TESTS FOR .validateConsensusConcordantsInput
# ==============================================================================

test_that(".validateConsensusConcordantsInput validates paired analysis requirements", {
    testData <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed)

    # Paired analysis requires exactly 2 dataframes
    expect_error(
        .validateConsensusConcordantsInput(list(testData), TRUE, 0.3, NULL),
        "Paired analysis requires two data frames"
    )

    # Paired analysis with 3 dataframes should error
    expect_error(
        .validateConsensusConcordantsInput(list(testData, testData, testData), TRUE, 0.3, NULL),
        "Paired analysis requires two data frames"
    )

    # Valid paired analysis with 2 dataframes
    expect_silent(
        .validateConsensusConcordantsInput(list(testData, testData), TRUE, 0.3, NULL)
    )
})

test_that(".validateConsensusConcordantsInput validates unpaired analysis requirements", {
    testData <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed)

    # Unpaired analysis requires exactly 1 dataframe
    expect_error(
        .validateConsensusConcordantsInput(list(testData, testData), FALSE, 0.3, NULL),
        "Unpaired analysis requires only one dataframe"
    )

    # Unpaired analysis with no dataframes should error
    expect_error(
        .validateConsensusConcordantsInput(list(), FALSE, 0.3, NULL),
        "Unpaired analysis requires only one dataframe"
    )

    # Valid unpaired analysis with 1 dataframe
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0.3, NULL)
    )
})

test_that(".validateConsensusConcordantsInput validates cutoff parameter", {
    testData <- createTestConcordants("CP", seed = .testSeed)

    # Non-numeric cutoff should error
    expect_error(
        .validateConsensusConcordantsInput(list(testData), FALSE, "0.3", NULL),
        "cutoff must be a single numeric value"
    )

    # Multiple cutoff values should error
    expect_error(
        .validateConsensusConcordantsInput(list(testData), FALSE, c(0.3, 0.5), NULL),
        "cutoff must be a single numeric value"
    )

    # Cutoff below 0 should error
    expect_error(
        .validateConsensusConcordantsInput(list(testData), FALSE, -0.1, NULL),
        "cutoff must be between 0 and 1"
    )

    # Cutoff above 1 should error
    expect_error(
        .validateConsensusConcordantsInput(list(testData), FALSE, 1.5, NULL),
        "cutoff must be between 0 and 1"
    )

    # Valid cutoff values should work
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0.321, NULL)
    )
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0L, NULL)
    )
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 1L, NULL)
    )
})

test_that(".validateConsensusConcordantsInput validates cellLine parameter", {
    testData <- createTestConcordants("CP", seed = .testSeed)

    # Non-character cellLine should error
    expect_error(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0.3, 123L),
        "cellLine must be a character vector or NULL"
    )

    # NULL cellLine should work
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0.3, NULL)
    )

    # Character vector cellLine should work
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0.3, "A375")
    )
    expect_silent(
        .validateConsensusConcordantsInput(list(testData), FALSE, 0.3, c("A375", "PC3"))
    )
})

test_that(".validateConsensusConcordantsInput validates dataframe content", {
    # Empty dataframe should error
    emptyData <- createTestConcordants("CP", nEntries = 0L)
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
    testData <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed)

    result <- .combineConcordantsData(list(testData))

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 5L)
    expect_identical(ncol(result), ncol(testData))
})

test_that(".combineConcordantsData combines multiple dataframes correctly", {
    testData1 <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed)
    testData2 <- createTestConcordants("CP", nEntries = 3L, seed = .testSeed + 1L)

    result <- .combineConcordantsData(list(testData1, testData2))

    expect_s3_class(result, "tbl_df")
    expect_identical(nrow(result), 8L) # 5 rows plus 3 rows
    expect_true(all(names(testData1) %in% names(result)))
})

# ==============================================================================
# TESTS FOR .filterByCellLine
# ==============================================================================

test_that(".filterByCellLine returns original data when cellLine is NULL", {
    testData <- createTestConcordants("CP", nEntries = 10L, seed = .testSeed)

    result <- .filterByCellLine(testData, NULL)

    expect_identical(result, testData)
})

test_that(".filterByCellLine filters by single cell line correctly", {
    testData <- createTestConcordants("CP",
        cellLines = c("A375", "PC3", "MCF7"),
        nEntries = 15L,
        seed = .testSeed
    )

    result <- .filterByCellLine(testData, "A375")

    expect_true(all(result[["cellline"]] == "A375"))
    expect_lt(nrow(result), nrow(testData))
})

test_that(".filterByCellLine filters by multiple cell lines correctly", {
    testData <- createTestConcordants("CP",
        cellLines = c("A375", "PC3", "MCF7", "HeLa"),
        nEntries = 20L,
        seed = .testSeed
    )

    result <- .filterByCellLine(testData, c("A375", "PC3"))

    expect_true(all(result[["cellline"]] %in% c("A375", "PC3")))
    expect_false(any(result[["cellline"]] %in% c("MCF7", "HeLa")))
})

test_that(".filterByCellLine returns empty result for non-existent cell lines", {
    testData <- createTestConcordants("CP", seed = .testSeed)

    result <- .filterByCellLine(testData, "NONEXISTENT")

    expect_identical(nrow(result), 0L)
    expect_s3_class(result, "tbl_df")
})

# ==============================================================================
# TESTS FOR .applySimilarityCutoff
# ==============================================================================

test_that(".applySimilarityCutoff filters by absolute similarity correctly", {
    testData <- tibble::tibble(
        similarity = c(0.5, -0.8, 0.2, -0.1, 0.6, -0.3),
        compound = paste0("COMP", 1L:6L)
    )

    result <- .applySimilarityCutoff(testData, 0.3)

    # Should keep |similarity| >= 0.3
    expect_identical(nrow(result), 4L) # 0.5, -0.8, 0.6, -0.3
    expect_true(all(abs(result[["similarity"]]) >= 0.3))
})

test_that(".applySimilarityCutoff handles zero cutoff", {
    testData <- tibble::tibble(
        similarity = c(0.5, -0.8, 0.0, -0.1),
        compound = paste0("COMP", 1L:4L)
    )

    result <- .applySimilarityCutoff(testData, 0L)

    # Should keep all non-negative absolute values
    expect_identical(nrow(result), nrow(testData))
})

test_that(".applySimilarityCutoff handles high cutoff", {
    testData <- tibble::tibble(
        similarity = c(0.5, -0.8, 0.2, -0.1),
        compound = paste0("COMP", 1L:4L)
    )

    result <- .applySimilarityCutoff(testData, 0.9)

    # Should return empty result
    expect_identical(nrow(result), 0L)
})

# ==============================================================================
# TESTS FOR .groupByTargetAndSelectMax
# ==============================================================================

test_that(".groupByTargetAndSelectMax selects maximum similarity per target", {
    testData <- tibble::tibble(
        compound = c("A", "A", "B", "B", "C"),
        similarity = c(0.5, 0.8, -0.3, -0.7, 0.4),
        cellline = c("A375", "PC3", "A375", "PC3", "MCF7"),
        signatureid = paste0("SIG", 1L:5L),
        treatment = c("A", "A", "B", "B", "C"),
        time = "24h",
        concentration = "10uM",
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .groupByTargetAndSelectMax(testData)

    # Should have one entry per compound with max |similarity|
    expect_identical(nrow(result), 3L)

    # Check that we got the maximum similarities
    expect_true(0.8 %in% result[["similarity"]]) # For compound A
    expect_true(-0.7 %in% result[["similarity"]]) # For compound B
    expect_true(0.4 %in% result[["similarity"]]) # For compound C
})

test_that(".groupByTargetAndSelectMax handles ties correctly", {
    testData <- tibble::tibble(
        compound = c("A", "A", "A"),
        similarity = c(0.5, -0.5, 0.3), # Two with |similarity| = 0.5
        cellline = c("A375", "PC3", "MCF7"),
        signatureid = paste0("SIG", 1L:3L),
        treatment = c("A", "A", "A"),
        time = "24h",
        concentration = "10uM",
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .groupByTargetAndSelectMax(testData)

    # Should keep both tied entries
    expect_identical(nrow(result), 2L)
    expect_true(all(abs(result$similarity) == 0.5))
})

test_that(".groupByTargetAndSelectMax works with treatment column", {
    testData <- tibble::tibble(
        treatment = c("GENE1", "GENE1", "GENE2"),
        similarity = c(0.5, 0.8, -0.3),
        cellline = c("A375", "PC3", "A375"),
        signatureid = paste0("SIG", 1L:3L),
        time = "24h",
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .groupByTargetAndSelectMax(testData)

    # Should have one entry per treatment with max |similarity|
    expect_identical(nrow(result), 2L)
    expect_true(0.8 %in% result$similarity)
    expect_true(-0.3 %in% result$similarity)
})

# ==============================================================================
# TESTS FOR .selectAndOrderResults
# ==============================================================================

test_that(".selectAndOrderResults selects correct columns", {
    testData <- tibble::tibble(
        signatureid = "SIG1",
        compound = "A",
        cellline = "A375",
        time = "24h",
        concentration = "10uM",
        similarity = 0.8,
        sig_direction = "Up",
        pValue = 0.01,
        extra_column = "EXTRA"
    )

    result <- .selectAndOrderResults(testData)

    expectedCols <- c(
        "signatureid", "compound", "cellline", "time",
        "concentration", "similarity", "sig_direction", "pValue"
    )
    expect_true(all(expectedCols %in% names(result)))
    expect_false("extra_column" %in% names(result))
})

test_that(".selectAndOrderResults orders by descending absolute similarity", {
    testData <- tibble::tibble(
        signatureid = paste0("SIG", 1L:5L),
        compound = paste0("COMP", 1L:5L),
        cellline = "A375",
        time = "24h",
        concentration = "10uM",
        similarity = c(0.3, -0.8, 0.5, -0.2, 0.7),
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .selectAndOrderResults(testData)

    # Should be ordered by descending |similarity|
    expectedOrder <- c(-0.8, 0.7, 0.5, 0.3, -0.2)
    expect_identical(result[["similarity"]], expectedOrder)
})

test_that(".selectAndOrderResults handles missing optional columns", {
    # Test with OE data (no concentration)
    testData <- tibble::tibble(
        signatureid = "SIG1",
        compound = "GENE1",
        cellline = "A375",
        time = "24h",
        similarity = 0.8,
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .selectAndOrderResults(testData)

    expect_s3_class(result, "tbl_df")
    expect_true("signatureid" %in% names(result))
    expect_true("similarity" %in% names(result))
})

# ==============================================================================
# TESTS FOR .applyTargetRenaming
# ==============================================================================

test_that(".applyTargetRenaming renames columns correctly for CP/KD libraries", {
    testData <- tibble::tibble(
        signatureid = "SIG1",
        compound = "A",
        cellline = "A375",
        time = "24h",
        concentration = "10uM",
        similarity = 0.8,
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .applyTargetRenaming(testData)

    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "Similarity",
        "SignatureDirection", "pValue"
    )
    expect_named(result, expectedCols)
})

test_that(".applyTargetRenaming renames columns correctly for OE library", {
    testData <- tibble::tibble(
        signatureid = "SIG1",
        treatment = "GENE1",
        cellline = "A375",
        time = "24h",
        similarity = 0.8,
        sig_direction = "Up",
        pValue = 0.01
    )

    result <- .applyTargetRenaming(testData)

    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "Similarity", "SignatureDirection", "pValue"
    )
    expect_named(result, expectedCols)
})

# ==============================================================================
# TESTS FOR .processConsensusPipeline
# ==============================================================================

test_that(".processConsensusPipeline processes data through complete pipeline", {
    testData <- createTestConcordants("CP", nEntries = 20L, seed = .testSeed) %>%
        dplyr::mutate(
            similarity = seq(-0.8, 0.8, length.out = 20L),
            cellline = rep(c("A375", "PC3", "MCF7"), length.out = 20L)
        )

    result <- .processConsensusPipeline(testData, 0.3, "A375")

    # Should be filtered by cell line
    expect_true(all(result[["TargetCellLine"]] == "A375"))

    # Should be filtered by cutoff
    expect_true(all(abs(result[["Similarity"]]) >= 0.3))

    # Should have renamed columns
    expect_true("TargetSignature" %in% names(result))
    expect_true("Target" %in% names(result))

    # Should be ordered by similarity
    if (nrow(result) > 1L) {
        similarities <- abs(result[["Similarity"]])
        expect_true(all(similarities[-length(similarities)] >= similarities[-1L]))
    }
})

test_that(".processConsensusPipeline handles NULL cellLine", {
    testData <- createTestConcordants("CP", nEntries = 10L, seed = .testSeed)

    result <- .processConsensusPipeline(testData, 0.3, NULL)

    expect_s3_class(result, "tbl_df")
    expect_true(all(abs(result$Similarity) >= 0.3))
})

# ==============================================================================
# TESTS FOR MAIN consensusConcordants FUNCTION
# ==============================================================================

test_that("consensusConcordants validates input correctly", {
    testData <- createTestConcordants("CP", seed = .testSeed)

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
    testData <- createTestConcordants("CP", nEntries = 20L, seed = .testSeed) %>%
        dplyr::mutate(similarity = seq(-0.8, 0.8, length.out = 20L))

    result <- consensusConcordants(testData, cutoff = 0.321)

    expect_s3_class(result, "tbl_df")
    expect_true(all(abs(result$Similarity) >= 0.321))
    expect_named(result, consensusConcordantsColumns("CP"))
})

test_that("consensusConcordants performs paired analysis correctly", {
    pairedData <- createPairedConcordants("CP", seed = .testSeed)

    result <- consensusConcordants(pairedData[[1L]], pairedData[[2L]], paired = TRUE)

    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
    expect_named(result, consensusConcordantsColumns("CP"))
})

test_that("consensusConcordants handles cell line filtering", {
    testData <- createTestConcordants("CP",
        cellLines = c("A375", "PC3", "MCF7"),
        nEntries = 30L,
        seed = .testSeed
    )
    # Ensure some entries for A375
    testData$cellline[1L:10L] <- "A375"

    result <- consensusConcordants(testData, cellLine = "A375")

    expect_true(all(result$TargetCellLine == "A375"))
})

test_that("consensusConcordants handles different cutoff values", {
    testData <- createTestConcordants("CP", nEntries = 20L, seed = .testSeed) %>%
        dplyr::mutate(similarity = seq(-0.8, 0.8, length.out = 20L))

    resultLow <- consensusConcordants(testData, cutoff = 0.1)
    resultHigh <- consensusConcordants(testData, cutoff = 0.7)

    expect_lt(nrow(resultHigh), nrow(resultLow))
    expect_true(all(abs(resultHigh$Similarity) >= 0.7))
})

test_that("consensusConcordants works with different library types", {
    # Test CP library
    cpData <- createTestConcordants("CP", seed = .testSeed)
    cpResult <- consensusConcordants(cpData)
    expect_true("TargetConcentration" %in% names(cpResult))

    # Test KD library
    kdData <- createTestConcordants("KD", seed = .testSeed)
    kdResult <- consensusConcordants(kdData)
    expect_false("TargetConcentration" %in% names(kdResult))

    # Test OE library
    oeData <- createTestConcordants("OE", seed = .testSeed)
    oeResult <- consensusConcordants(oeData, cutoff = 0.2)
    expect_false("TargetConcentration" %in% names(oeResult))
})

test_that("consensusConcordants handles empty results gracefully", {
    testData <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed) %>%
        dplyr::mutate(similarity = c(0.1, -0.1, 0.05, -0.05, 0.0))

    result <- consensusConcordants(testData, cutoff = 0.9)

    expect_equal(nrow(result), 0L)
    expect_s3_class(result, "tbl_df")
})

test_that("consensusConcordants maintains data integrity", {
    testData <- createTestConcordants("CP", seed = .testSeed)

    result <- consensusConcordants(testData)

    # Check data types
    expect_true(is.character(result$TargetSignature))
    expect_true(is.character(result$Target))
    expect_true(is.numeric(result$Similarity))
    expect_true(is.numeric(result$pValue))

    # Check no missing critical values
    expect_false(any(is.na(result$Target)))
    expect_false(any(is.infinite(result$Similarity)))
})
# ==============================================================================

test_that("consensusConcordants validates input arguments correctly", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed)

    # Test paired analysis with wrong number of arguments
    expect_error(
        consensusConcordants(testConcordants, paired = TRUE),
        "Paired analysis requires two data frames"
    )

    # Test unpaired analysis with too many arguments
    expect_error(
        consensusConcordants(testConcordants, testConcordants, paired = FALSE),
        "Unpaired analysis requires only one dataframe"
    )

    # Valid single dataframe should work
    expect_silent({
        result <- consensusConcordants(testConcordants, paired = FALSE)
    })

    # Valid paired dataframes should work
    pairedData <- createPairedConcordants("CP", seed = .testSeed)
    expect_silent({
        result <- consensusConcordants(pairedData[[1L]], pairedData[[2L]], paired = TRUE)
    })
})

test_that("consensusConcordants handles empty input correctly", {
    emptyData <- createTestConcordants("CP", nEntries = 0L)

    expect_error(consensusConcordants(emptyData, paired = FALSE), "All input dataframes must be non-empty")
})

test_that("consensusConcordants handles NULL cellLine parameter", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed)

    # Should work without cellLine filtering
    expect_silent({
        result <- consensusConcordants(testConcordants, cellLine = NULL)
    })

    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
})

# ==============================================================================
# TESTS FOR CUTOFF FILTERING
# ==============================================================================

test_that("consensusConcordants properly handles similarity threshold", {
    testConcordants <- createTestConcordants("CP", nEntries = 20L, seed = .testSeed) %>%
        dplyr::mutate(
            # Ensure we have data above and below cutoff
            similarity = c(rep(0.5, 5L), rep(0.2, 5L), rep(-0.5, 5L), rep(-0.2, 5L))
        )

    result <- consensusConcordants(testConcordants, cutoff = 0.321)

    # All results should meet cutoff
    expect_true(all(abs(result[["Similarity"]]) >= 0.321))

    # Should filter out low similarity values
    expect_lt(nrow(result), nrow(testConcordants))

    # Check column structure
    expect_named(result, consensusConcordantsColumns("CP"))
})

test_that("consensusConcordants respects different cutoff values", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed) %>%
        dplyr::mutate(similarity = seq(-0.8, 0.8, length.out = nrow(.)))

    # Test with low cutoff
    resultLow <- consensusConcordants(testConcordants, cutoff = 0.1)

    # Test with high cutoff
    resultHigh <- consensusConcordants(testConcordants, cutoff = 0.7)

    # High cutoff should return fewer results
    expect_lt(nrow(resultHigh), nrow(resultLow))
    expect_true(all(abs(resultHigh[["Similarity"]]) >= 0.7))
})

# ==============================================================================
# TESTS FOR CELL LINE FILTERING
# ==============================================================================

test_that("consensusConcordants properly handles single cell line filtering", {
    testConcordants <- createTestConcordants("CP",
        cellLines = c("A375", "PC3", "MCF7"),
        nEntries = 15L,
        seed = .testSeed
    )

    result <- consensusConcordants(testConcordants, cellLine = "A375")

    # All results should be from A375
    expect_true(all(result[["TargetCellLine"]] == "A375"))

    # Should have fewer results than total
    expect_lt(nrow(result), nrow(testConcordants))

    # Check column structure
    expect_named(result, consensusConcordantsColumns("CP"))
})

test_that("consensusConcordants properly handles multiple cell line filtering", {
    testConcordants <- createTestConcordants("CP",
        cellLines = c("A375", "PC3", "MCF7", "HeLa"),
        nEntries = 20L,
        seed = .testSeed
    )

    result <- consensusConcordants(testConcordants, cellLine = c("A375", "PC3"))

    # All results should be from specified cell lines
    expect_true(all(result[["TargetCellLine"]] %in% c("A375", "PC3")))

    # Should exclude other cell lines
    expect_false(any(result[["TargetCellLine"]] %in% c("MCF7", "HeLa")))

    # Check column structure
    expect_named(result, consensusConcordantsColumns("CP"))
})

test_that("consensusConcordants handles non-existent cell lines gracefully", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed)

    result <- consensusConcordants(testConcordants, cellLine = "NONEXISTENT")

    # Should return empty result
    expect_equal(nrow(result), 0L)
    expect_s3_class(result, "tbl_df")
})

# ==============================================================================
# TESTS FOR PAIRED ANALYSIS
# ==============================================================================

test_that("consensusConcordants properly handles paired analysis", {
    pairedData <- createPairedConcordants("CP", seed = .testSeed)
    upConcordants <- pairedData[[1L]]
    downConcordants <- pairedData[[2L]]

    result <- consensusConcordants(upConcordants, downConcordants, paired = TRUE)

    # Check structure
    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
    expect_named(result, consensusConcordantsColumns("CP"))

    # Should combine data from both inputs
    expect_true(any(result[["SignatureDirection"]] == "Up"))
    expect_true(any(result[["SignatureDirection"]] == "Down"))

    # Should not have duplicate targets (due to grouping)
    duplicatedTargets <- duplicated(result[["Target"]])
    expect_false(all(duplicatedTargets))
})

test_that("consensusConcordants paired analysis with cell line filtering works", {
    pairedData <- createPairedConcordants("CP", seed = .testSeed)
    # Ensure some data has A375 cell line
    pairedData[[1L]]$cellline[1L:2L] <- "A375"
    pairedData[[2L]]$cellline[1L:2L] <- "A375"

    result <- consensusConcordants(
        pairedData[[1L]], pairedData[[2L]],
        paired = TRUE,
        cellLine = "A375"
    )

    # All results should be from A375
    expect_true(all(result[["TargetCellLine"]] == "A375"))

    # Check structure
    expect_named(result, consensusConcordantsColumns("CP"))
    expect_s3_class(result, "tbl_df")
})

# ==============================================================================
# TESTS FOR DATA PROCESSING AND TRANSFORMATION
# ==============================================================================

test_that("consensusConcordants correctly groups and filters by maximum similarity", {
    # Create test data with duplicate compounds but different similarities
    testConcordants <- tibble::tibble(
        signatureid = c("SIG1", "SIG2", "SIG3", "SIG4"),
        compound = c("COMP1", "COMP1", "COMP2", "COMP2"),
        cellline = c("A375", "A375", "PC3", "PC3"),
        time = c("24h", "24h", "24h", "24h"),
        concentration = c("10uM", "10uM", "10uM", "10uM"),
        similarity = c(0.5, 0.8, -0.3, -0.7), # COMP1 max: 0.8, COMP2 max: -0.7
        sig_direction = c("Up", "Up", "Down", "Down"),
        pValue = c(0.01, 0.02, 0.03, 0.04),
        sig_type = rep("Chemical Perturbagen", 4L),
        nGenes = rep(978L, 4L),
        lincsPertID = paste0("LSM-", 1L:4L),
        GeneTargets = rep("NA", 4L),
        `_row` = paste0("SIG", 1L:4L)
    )

    result <- consensusConcordants(testConcordants, cutoff = 0.2)

    # Should only keep maximum similarity for each compound
    expect_equal(nrow(result), 2L) # One for each compound

    # Check that we got the maximum similarities
    comp1_result <- result[result$Target == "COMP1", ]
    comp2_result <- result[result$Target == "COMP2", ]

    expect_equal(comp1_result$Similarity, 0.8)
    expect_equal(comp2_result$Similarity, -0.7)
})

test_that("consensusConcordants correctly renames columns using targetRename", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed)

    result <- consensusConcordants(testConcordants)

    # Check that columns are renamed correctly
    expectedCols <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "Similarity",
        "SignatureDirection", "pValue"
    )

    expect_named(result, expectedCols)
})

test_that("consensusConcordants correctly sorts by similarity magnitude", {
    testConcordants <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed) %>%
        dplyr::mutate(
            similarity = c(0.3, -0.8, 0.5, -0.2, 0.7),
        )

    result <- consensusConcordants(testConcordants, cutoff = 0.1)

    # Should be sorted by descending absolute similarity
    expectedOrder <- c(-0.8, 0.7, 0.5, 0.3, -0.2)
    expect_equal(result$Similarity, expectedOrder)
})

# ==============================================================================
# TESTS FOR DIFFERENT LIBRARY TYPES
# ==============================================================================

test_that("consensusConcordants handles CP library correctly", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed)

    result <- consensusConcordants(testConcordants)

    # CP should have concentration column
    expect_true("TargetConcentration" %in% names(result))
    expect_named(result, consensusConcordantsColumns("CP"))
})

test_that("consensusConcordants handles KD library correctly", {
    testConcordants <- createTestConcordants("KD", seed = .testSeed)

    result <- consensusConcordants(testConcordants)

    # KD should have concentration column (renamed from treatment)
    expect_false("TargetConcentration" %in% names(result))
    expect_named(result, consensusConcordantsColumns("KD"))
})

test_that("consensusConcordants handles OE library correctly", {
    testConcordants <- createTestConcordants("OE", seed = .testSeed)

    result <- consensusConcordants(testConcordants, cutoff = 0.4)

    # OE should NOT have concentration column
    expect_false("TargetConcentration" %in% names(result))
    expect_named(result, consensusConcordantsColumns("OE"))
    expect_s3_class(result, "tbl_df")
})

# ==============================================================================
# TESTS FOR EDGE CASES AND ERROR CONDITIONS
# ==============================================================================

test_that("consensusConcordants handles missing required columns gracefully", {
    # Test with missing similarity column
    incompleteData <- tibble::tibble(
        signatureid = "SIG1",
        compound = "COMP1",
        cellline = "A375"
    )

    expect_error({
        consensusConcordants(incompleteData)
    })
})

test_that("consensusConcordants handles NA values appropriately", {
    testConcordants <- createTestConcordants("CP", nEntries = 5L, seed = .testSeed) %>%
        dplyr::mutate(
            similarity = c(0.5, NA, 0.3, -0.4, 0.6),
            pValue = c(0.01, 0.02, NA, 0.04, 0.05)
        )

    result <- consensusConcordants(testConcordants, cutoff = 0.2)

    # Should handle NA values without error
    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
})

test_that("consensusConcordants maintains data integrity", {
    testConcordants <- createTestConcordants("CP", seed = .testSeed)

    result <- consensusConcordants(testConcordants)

    # Check data types
    expect_true(is.character(result$TargetSignature))
    expect_true(is.character(result$Target))
    expect_true(is.character(result$TargetCellLine))
    expect_true(is.numeric(result$Similarity))
    expect_true(is.numeric(result$pValue))

    # Check no infinite or missing critical values
    expect_false(any(is.infinite(result$Similarity)))
    expect_false(any(is.na(result$Target)))
})

# ==============================================================================
# INTEGRATION TESTS USING FIXTURES
# ==============================================================================

test_that("consensusConcordants works with test fixtures", {
    # Test with CP fixture
    testConcordants <- getTestFixture("concordants", "CP", seed = .testSeed)

    result <- consensusConcordants(testConcordants, cutoff = 0.321)

    expect_s3_class(result, "tbl_df")
    expect_true(all(abs(result[["Similarity"]]) >= 0.321))
    expect_named(result, consensusConcordantsColumns("CP"))
})
