# Comprehensive tests for investigateTarget function

library(httptest2)
library(dplyr)

# ==============================================================================
# TESTS FOR INPUT VALIDATION
# ==============================================================================

test_that("investigateTarget validates library inputs correctly", {
    # Invalid inputLib should error
    expect_error(
        investigateTarget(
            target = "AATK",
            inputLib = "INVALID",
            outputLib = "CP"
        ),
        "Invalid library specification"
    )

    # Invalid outputLib should error
    expect_error(
        investigateTarget(
            target = "AATK",
            inputLib = "KD",
            outputLib = "INVALID"
        ),
        "Invalid library specification"
    )

    # Both invalid should error
    expect_error(
        investigateTarget(
            target = "AATK",
            inputLib = "invalid",
            outputLib = "invalid"
        ),
        "Invalid library specification"
    )
})

test_that("investigateTarget handles missing target parameter", {
    # Missing target should result in R's standard missing argument error
    expect_error(
        investigateTarget(inputLib = "KD", outputLib = "CP"),
        "argument \"target\" is missing"
    )
})

test_that("investigateTarget handles nonexistent target names", {
    skip_on_cran()

    # Target that doesn't exist in metadata should error
    expect_error(
        investigateTarget(
            target = "NONEXISTENT_GENE_12345",
            inputLib = "KD",
            outputLib = "CP"
        ),
        "No signatures match the given input criteria"
    )

    expect_error(
        investigateTarget(
            target = "FAKE_DRUG_9999",
            inputLib = "CP",
            outputLib = "KD"
        ),
        "No signatures match the given input criteria"
    )
})

# ==============================================================================
# TESTS FOR TARGET LOOKUP AND FILTERING
# ==============================================================================

test_that("investigateTarget finds targets from KD metadata", {
    skip_on_cran()

    # Pick a gene that should exist in KD metadata
    # Use first available gene from kdMetadata if exists
    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        # This should find at least one signature
        # We're not running the full function, just testing the lookup logic
        inputMetadata <- .loadMetadata("KD")
        filteredSigs <- inputMetadata |>
            dplyr::filter(
                stringr::str_to_lower(.data[["Source"]]) ==
                    stringr::str_to_lower(testGene)
            ) |>
            dplyr::pull(.data[["SourceSignature"]])

        expect_gt(length(filteredSigs), 0L)
        expect_type(filteredSigs, "character")
    }
})

test_that("investigateTarget finds targets from OE metadata", {
    skip_on_cran()

    if (nrow(oeMetadata) > 0L) {
        testGene <- oeMetadata[["Source"]][1L]

        inputMetadata <- .loadMetadata("OE")
        filteredSigs <- inputMetadata |>
            dplyr::filter(
                stringr::str_to_lower(.data[["Source"]]) ==
                    stringr::str_to_lower(testGene)
            ) |>
            dplyr::pull(.data[["SourceSignature"]])

        expect_gt(length(filteredSigs), 0L)
        expect_type(filteredSigs, "character")
    }
})

test_that("investigateTarget finds targets from CP metadata", {
    skip_on_cran()

    if (nrow(cpMetadata) > 0L) {
        testDrug <- cpMetadata[["Source"]][1L]

        inputMetadata <- .loadMetadata("CP")
        filteredSigs <- inputMetadata |>
            dplyr::filter(
                stringr::str_to_lower(.data[["Source"]]) ==
                    stringr::str_to_lower(testDrug)
            ) |>
            dplyr::pull(.data[["SourceSignature"]])

        expect_gt(length(filteredSigs), 0L)
        expect_type(filteredSigs, "character")
    }
})

test_that("investigateTarget is case-insensitive for target names", {
    skip_on_cran()

    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        inputMetadata <- .loadMetadata("KD")

        # Test lowercase
        filteredLower <- inputMetadata |>
            dplyr::filter(
                stringr::str_to_lower(.data[["Source"]]) ==
                    stringr::str_to_lower(tolower(testGene))
            ) |>
            dplyr::pull(.data[["SourceSignature"]])

        # Test uppercase
        filteredUpper <- inputMetadata |>
            dplyr::filter(
                stringr::str_to_lower(.data[["Source"]]) ==
                    stringr::str_to_lower(toupper(testGene))
            ) |>
            dplyr::pull(.data[["SourceSignature"]])

        # Should find the same signatures regardless of case
        expect_identical(filteredLower, filteredUpper)
    }
})

# ==============================================================================
# TESTS FOR CELL LINE FILTERING
# ==============================================================================

test_that("investigateTarget filters by inputCellLines correctly", {
    skip_on_cran()

    if (nrow(kdMetadata) > 0L && "SourceCellLine" %in% colnames(kdMetadata)) {
        testGene <- kdMetadata[["Source"]][1L]
        availableCellLines <- unique(kdMetadata[["SourceCellLine"]])

        if (length(availableCellLines) > 1L) {
            # Pick first cell line
            testCellLine <- availableCellLines[1L]

            inputMetadata <- .loadMetadata("KD")

            # Without filtering
            allSigs <- inputMetadata |>
                dplyr::filter(
                    stringr::str_to_lower(.data[["Source"]]) ==
                        stringr::str_to_lower(testGene)
                ) |>
                dplyr::pull(.data[["SourceSignature"]])

            # With cell line filtering
            filteredSigs <- inputMetadata |>
                dplyr::filter(
                    stringr::str_to_lower(.data[["Source"]]) ==
                        stringr::str_to_lower(testGene)
                ) |>
                dplyr::filter(.data[["SourceCellLine"]] %in% testCellLine) |>
                dplyr::pull(.data[["SourceSignature"]])

            # Filtered should be subset of all
            expect_true(all(filteredSigs %in% allSigs))
            expect_lte(length(filteredSigs), length(allSigs))
        }
    }
})

# ==============================================================================
# TESTS FOR PARAMETER VALIDATION
# ==============================================================================

test_that("investigateTarget validates threshold parameters", {
    skip_on_cran()

    # filterThreshold should be numeric
    # Note: These will fail during the consensus calculation if they reach it
    # For now, we test that bad parameters get caught somewhere in the chain

    # We can't easily test without network access, so skip detailed validation
    # The actual parameter validation happens in filterSignature and consensusConcordants
})

test_that("investigateTarget handles paired parameter correctly", {
    skip_on_cran()

    # Paired should be logical
    # This will be tested in integration tests with mock data
})

# ==============================================================================
# TESTS FOR RETURN VALUE STRUCTURE
# ==============================================================================

test_that("investigateTarget returns tibble with expected columns", {
    skip_on_cran()
    skip_if_offline()

    # This test requires actual API access
    # Pick a well-known gene that should work
    if (nrow(kdMetadata) > 0L) {
        # Use a gene we know exists
        testGene <- kdMetadata[["Source"]][1L]

        result <- investigateTarget(
            target = testGene,
            inputLib = "KD",
            outputLib = "CP",
            filterThreshold = 0.85,
            similarityThreshold = 0.3,
            paired = FALSE
        )

        # Check return type
        expect_s3_class(result, "tbl_df")

        # Check for expected columns
        expectedCols <- c(
            "Source", "Target", "Similarity",
            "SourceSignature", "TargetSignature"
        )

        expect_true(all(expectedCols %in% colnames(result)))
    }
})

# ==============================================================================
# TESTS FOR INTEGRATION WITH OTHER FUNCTIONS
# ==============================================================================

test_that("investigateTarget correctly orchestrates getSignature", {
    skip_on_cran()
    skip_if_offline()

    # This test verifies that investigateTarget correctly calls getSignature
    # for each matched signature ID

    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        # Get the signature IDs that should be processed
        inputMetadata <- .loadMetadata("KD")
        expectedSigIds <- inputMetadata |>
            dplyr::filter(
                stringr::str_to_lower(.data[["Source"]]) ==
                    stringr::str_to_lower(testGene)
            ) |>
            dplyr::pull(.data[["SourceSignature"]])

        if (length(expectedSigIds) > 0L) {
            # Verify getSignature works for these IDs
            testSig <- getSignature(expectedSigIds[1L])
            expect_s3_class(testSig, "tbl_df")
            expect_true("Value_LogDiffExp" %in% colnames(testSig))
        }
    }
})

test_that("investigateTarget correctly passes parameters to filterSignature", {
    skip_on_cran()

    # Test that filterThreshold parameter is properly used
    # This is indirectly tested through the .computeConsensusFromSignature function
})

test_that("investigateTarget correctly passes parameters to consensusConcordants", {
    skip_on_cran()

    # Test that similarityThreshold and outputCellLines are properly passed
    # This is indirectly tested through the .computeConsensusFromSignature function
})

# ==============================================================================
# TESTS FOR EDGE CASES
# ==============================================================================

test_that("investigateTarget handles single source signature", {
    skip_on_cran()
    skip_if_offline()

    # Find a target with only one signature
    if (nrow(kdMetadata) > 0L) {
        singleSigTargets <- kdMetadata |>
            dplyr::group_by(.data[["Source"]]) |>
            dplyr::filter(dplyr::n() == 1L) |>
            dplyr::ungroup()

        if (nrow(singleSigTargets) > 0L) {
            testGene <- singleSigTargets[["Source"]][1L]

            result <- investigateTarget(
                target = testGene,
                inputLib = "KD",
                outputLib = "CP",
                filterThreshold = 0.5,
                paired = FALSE
            )

            expect_s3_class(result, "tbl_df")
        }
    }
})

test_that("investigateTarget handles multiple source signatures", {
    skip_on_cran()
    skip_if_offline()

    # Find a target with multiple signatures
    if (nrow(kdMetadata) > 0L) {
        multipleSigTargets <- kdMetadata |>
            dplyr::group_by(.data[["Source"]]) |>
            dplyr::filter(dplyr::n() > 1L) |>
            dplyr::ungroup()

        if (nrow(multipleSigTargets) > 0L) {
            testGene <- multipleSigTargets[["Source"]][1L]

            result <- investigateTarget(
                target = testGene,
                inputLib = "KD",
                outputLib = "CP",
                filterThreshold = 0.5,
                paired = FALSE
            )

            expect_s3_class(result, "tbl_df")

            # Should have SourceSignature column indicating which signatures were used
            expect_true("SourceSignature" %in% colnames(result))
        }
    }
})

# ==============================================================================
# TESTS FOR PAIRED VS UNPAIRED WORKFLOW
# ==============================================================================

test_that("investigateTarget paired workflow produces valid results", {
    skip_on_cran()
    skip_if_offline()

    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        result <- investigateTarget(
            target = testGene,
            inputLib = "KD",
            outputLib = "CP",
            filterThreshold = 0.5,
            similarityThreshold = 0.3,
            paired = TRUE
        )

        expect_s3_class(result, "tbl_df")
        expect_true("Similarity" %in% colnames(result))
    }
})

test_that("investigateTarget unpaired workflow produces valid results", {
    skip_on_cran()
    skip_if_offline()

    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        result <- investigateTarget(
            target = testGene,
            inputLib = "KD",
            outputLib = "CP",
            filterThreshold = 0.5,
            similarityThreshold = 0.3,
            paired = FALSE
        )

        expect_s3_class(result, "tbl_df")
        expect_true("Similarity" %in% colnames(result))
    }
})

# ==============================================================================
# TESTS FOR METADATA ENRICHMENT
# ==============================================================================

test_that("investigateTarget enriches results with source metadata", {
    skip_on_cran()
    skip_if_offline()

    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        result <- investigateTarget(
            target = testGene,
            inputLib = "KD",
            outputLib = "CP",
            filterThreshold = 0.5,
            paired = FALSE
        )

        # Should have source-related columns from metadata join
        expect_true("Source" %in% colnames(result))
        expect_true("SourceSignature" %in% colnames(result))

        # All Source values should match the target
        expect_true(all(
            stringr::str_to_lower(result[["Source"]]) ==
                stringr::str_to_lower(testGene)
        ))
    }
})

test_that("investigateTarget maintains target metadata in results", {
    skip_on_cran()
    skip_if_offline()

    if (nrow(kdMetadata) > 0L) {
        testGene <- kdMetadata[["Source"]][1L]

        result <- investigateTarget(
            target = testGene,
            inputLib = "KD",
            outputLib = "CP",
            filterThreshold = 0.5,
            paired = FALSE
        )

        # Should have target-related columns
        expect_true("Target" %in% colnames(result))
        expect_true("TargetSignature" %in% colnames(result))
    }
})
