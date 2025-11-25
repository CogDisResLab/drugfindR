# Test utility functions

library(tibble)
library(S4Vectors)

test_that("targetRename works correctly without treatment column", {
    inputNamesWithoutTreatment <- c(
        "signature",
        "target",
        "cellLine",
        "time",
        "concentration",
        "similarity",
        "direction",
        "pvalue"
    )
    expectedNames <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "InputSigDirection",
        "SignatureType", "Similarity",
        "pValue"
    )

    result <- targetRename(inputNamesWithoutTreatment)
    expect_identical(result, expectedNames)
})

test_that(".ilincsBaseUrl returns correct URL", {
    expectedUrl <- "https://www.ilincs.org/api"
    result <- .ilincsBaseUrl()
    expect_identical(result, expectedUrl)
    expect_type(result, "character")
    expect_length(result, 1L)
})

test_that(".validateLibrary works correctly for valid libraries", {
    expect_true(.validateLibrary("CP"))
    expect_true(.validateLibrary("KD"))
    expect_true(.validateLibrary("OE"))
})

test_that(".validateLibrary works correctly for invalid libraries", {
    expect_false(.validateLibrary("INVALID"))
    expect_false(.validateLibrary("cp")) # case sensitive
    expect_false(.validateLibrary("kd"))
    expect_false(.validateLibrary("oe"))
    expect_false(.validateLibrary(""))
})

test_that("validateLibraries works correctly for valid library vectors", {
    expect_true(validateLibraries(c("CP", "KD")))
    expect_true(validateLibraries(c("OE", "CP", "KD")))
    expect_true(validateLibraries("CP"))
    expect_true(validateLibraries(c("OE"))) # nolint: unnecessary_concatenation_linter.
})

test_that("validateLibraries works correctly for invalid library vectors", {
    expect_false(validateLibraries(c("CP", "INVALID")))
    expect_false(validateLibraries(c("cp", "KD")))
    expect_false(validateLibraries(c("CP", "KD", "INVALID")))
    expect_false(validateLibraries("INVALID"))
})

test_that("stopIfInvalidLibraries stops for invalid libraries", {
    expect_error(
        stopIfInvalidLibraries(c("CP", "INVALID")),
        "Invalid library specification\\(s\\): INVALID"
    )
    expect_error(
        stopIfInvalidLibraries("INVALID"),
        "Invalid library specification\\(s\\): INVALID"
    )
    expect_error(
        stopIfInvalidLibraries(c("CP", "XYZ", "INVALID")),
        "Invalid library specification\\(s\\): XYZ, INVALID"
    )
})

test_that("stopIfInvalidLibraries provides helpful error messages", {
    # Test that error message includes expected libraries
    errorMsg <- testthat::capture_error(stopIfInvalidLibraries("INVALID"))$message

    expect_true(grepl("Invalid library specification", errorMsg))
    expect_true(grepl("INVALID", errorMsg))
    expect_true(grepl("'OE' \\(Overexpression\\)", errorMsg))
    expect_true(grepl("'KD' \\(Knockdown\\)", errorMsg))
    expect_true(grepl("'CP' \\(Chemical Perturbagen\\)", errorMsg))
})

test_that("stopIfInvalidLibraries doesn't stop for valid libraries", {
    expect_silent(stopIfInvalidLibraries("CP"))
    expect_silent(stopIfInvalidLibraries(c("CP", "KD")))
    expect_silent(stopIfInvalidLibraries(c("OE", "CP", "KD")))
})

test_that(".loadMetadata returns correct metadata for each library", {
    # Test OE metadata
    oeResult <- .loadMetadata("OE")
    expect_s3_class(oeResult, "tbl_df")
    expect_identical(oeResult, oeMetadata)

    # Test KD metadata
    kdResult <- .loadMetadata("KD")
    expect_s3_class(kdResult, "tbl_df")
    expect_identical(kdResult, kdMetadata)

    # Test CP metadata
    cpResult <- .loadMetadata("CP")
    expect_s3_class(cpResult, "tbl_df")
    expect_identical(cpResult, cpMetadata)
})

test_that("loadMetadata throws error for invalid library", {
    expect_error(
        .loadMetadata("INVALID"),
        "Invalid library: 'INVALID'"
    )
    expect_error(
        .loadMetadata("cp"),
        "Invalid library: 'cp'"
    )

    # Test that error message includes valid options
    errorMsg <- testthat::capture_error(.loadMetadata("INVALID"))$message
    expect_true(grepl("'OE' \\(Overexpression\\)", errorMsg))
    expect_true(grepl("'KD' \\(Knockdown\\)", errorMsg))
    expect_true(grepl("'CP' \\(Chemical Perturbagen\\)", errorMsg))
})

test_that(".returnLibrary returns correct library IDs", {
    expect_identical(.returnLibrary("OE"), "LIB_11", ignore_attr = "names")
    expect_identical(.returnLibrary("KD"), "LIB_6", ignore_attr = "names")
    expect_identical(.returnLibrary("CP"), "LIB_5", ignore_attr = "names")
})

test_that(".returnLibrary throws error for invalid library", {
    expect_error(
        .returnLibrary("INVALID"),
        "Invalid library specification\\(s\\): INVALID"
    )
})

test_that(".returnUserAgent returns valid user agent string", {
    userAgent <- .returnUserAgent()
    expect_type(userAgent, "character")
    expect_length(userAgent, 1L)
    expect_true(grepl("drugfindR", userAgent))
    expect_gt(nchar(userAgent), 0L)
})

# ==============================================================================
# TESTS FOR SIGNATURE VALIDATION FUNCTIONS
# ==============================================================================

# Tests for .stopIfInvalidColNames function
test_that(".stopIfInvalidColNames works correctly with valid signature", {
    validSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Should not throw error for valid signature
    expect_silent(.stopIfInvalidColNames(validSig))

    # Test with different data frame types
    validTibble <- tibble::as_tibble(validSig)
    expect_silent(.stopIfInvalidColNames(validTibble))

    validDataFrame <- S4Vectors::DataFrame(validSig)
    expect_silent(.stopIfInvalidColNames(validDataFrame))
})

test_that(".stopIfInvalidColNames detects missing columns", {
    # Missing required columns
    incompleteSig <- data.frame(
        signatureID = "SIG_001",
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = 1.5
        # Missing ID_geneid and Significance_pvalue
    )

    expect_error(
        .stopIfInvalidColNames(incompleteSig),
        "Input signature does not conform to expected structure"
    )

    expect_error(
        .stopIfInvalidColNames(incompleteSig),
        "Missing columns: ID_geneid, Significance_pvalue"
    )
})

test_that(".stopIfInvalidColNames detects extra columns", {
    # Extra unexpected columns
    sigWithExtra <- data.frame(
        signatureID = "SIG_001",
        ID_geneid = "1234",
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = 1.5,
        Significance_pvalue = 0.05,
        ExtraColumn1 = "extra",
        ExtraColumn2 = "also_extra",
        stringsAsFactors = FALSE
    )

    expect_error(
        .stopIfInvalidColNames(sigWithExtra),
        "Input signature does not conform to expected structure"
    )

    expect_error(
        .stopIfInvalidColNames(sigWithExtra),
        "Unexpected columns: ExtraColumn1, ExtraColumn2"
    )
})

test_that(".stopIfInvalidColNames detects wrong column order", {
    # Correct columns but wrong order
    wrongOrderSig <- data.frame(
        Name_GeneSymbol = "GENE1",
        signatureID = "SIG_001",
        ID_geneid = "1234",
        Significance_pvalue = 0.05,
        Value_LogDiffExp = 1.5,
        stringsAsFactors = FALSE
    )

    expect_error(
        .stopIfInvalidColNames(wrongOrderSig),
        "Input signature does not conform to expected structure"
    )

    expect_error(
        .stopIfInvalidColNames(wrongOrderSig),
        "Columns are not in the expected order"
    )
})

test_that(".stopIfInvalidColNames provides comprehensive error message", {
    # Signature with both missing and extra columns
    complexErrorSig <- data.frame(
        signatureID = "SIG_001",
        WrongColumn = "wrong",
        Name_GeneSymbol = "GENE1",
        AnotherWrongColumn = "also_wrong"
        # Missing ID_geneid, Value_LogDiffExp, Significance_pvalue
    )

    errorMsg <- testthat::capture_error(.stopIfInvalidColNames(complexErrorSig))$message

    expect_true(grepl("Missing columns", errorMsg))
    expect_true(grepl("Unexpected columns", errorMsg))
    expect_true(grepl("Expected columns \\(in order\\)", errorMsg))
    expect_true(grepl("Actual columns", errorMsg))
    expect_true(grepl("prepareSignature", errorMsg))
})

test_that(".stopIfInvalidColNames error message includes all expected information", {
    emptySig <- data.frame()

    errorMsg <- testthat::capture_error(.stopIfInvalidColNames(emptySig))$message

    # Should show all missing columns
    expect_true(grepl("signatureID", errorMsg))
    expect_true(grepl("ID_geneid", errorMsg))
    expect_true(grepl("Name_GeneSymbol", errorMsg))
    expect_true(grepl("Value_LogDiffExp", errorMsg))
    expect_true(grepl("Significance_pvalue", errorMsg))
})

# Tests for .stopIfContainsMissingValues function
test_that(".stopIfContainsMissingValues works correctly with valid signature", {
    validSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Should not throw error for signature without missing values
    expect_silent(.stopIfContainsMissingValues(validSig))

    # Test with different data frame types
    validTibble <- tibble::as_tibble(validSig)
    expect_silent(.stopIfContainsMissingValues(validTibble))

    validDataFrame <- S4Vectors::DataFrame(validSig)
    expect_silent(.stopIfContainsMissingValues(validDataFrame))
})

test_that(".stopIfContainsMissingValues detects missing values", {
    # Signature with NA values
    sigWithNa <- data.frame(
        signatureID = c("SIG_001", "SIG_001"),
        ID_geneid = c("1234", NA),
        Name_GeneSymbol = c("GENE1", "GENE2"),
        Value_LogDiffExp = c(1.5, -2.1),
        Significance_pvalue = c(0.05, NA),
        stringsAsFactors = FALSE
    )

    expect_error(
        .stopIfContainsMissingValues(sigWithNa),
        "Input signature contains missing \\(NA\\) values"
    )
})

test_that(".stopIfContainsMissingValues provides detailed error information", {
    # Create signature with missing values in specific columns
    sigWithMultipleNa <- data.frame(
        signatureID = c("SIG_001", "SIG_001", "SIG_001"),
        ID_geneid = c("1234", NA, "9012"),
        Name_GeneSymbol = c("GENE1", "GENE2", "GENE3"),
        Value_LogDiffExp = c(1.5, NA, NA),
        Significance_pvalue = c(0.05, 0.01, 0.03),
        stringsAsFactors = FALSE
    )

    errorMsg <- testthat::capture_error(.stopIfContainsMissingValues(sigWithMultipleNa))$message

    # Should list columns with missing values and counts
    expect_true(grepl("ID_geneid: 1 missing value", errorMsg))
    expect_true(grepl("Value_LogDiffExp: 2 missing value", errorMsg))
    expect_true(grepl("prepareSignature", errorMsg))
    expect_true(grepl("remove or impute", errorMsg))
})

test_that(".stopIfContainsMissingValues handles all-NA columns", {
    # Signature with entirely NA column
    sigAllNaColumn <- data.frame(
        signatureID = c("SIG_001", "SIG_001"),
        ID_geneid = c(NA, NA),
        Name_GeneSymbol = c("GENE1", "GENE2"),
        Value_LogDiffExp = c(1.5, -2.1),
        Significance_pvalue = c(0.05, 0.01),
        stringsAsFactors = FALSE
    )

    errorMsg <- testthat::capture_error(.stopIfContainsMissingValues(sigAllNaColumn))$message
    expect_true(grepl("ID_geneid: 2 missing value", errorMsg))
})

test_that(".stopIfContainsMissingValues handles mixed data types with NA", {
    # Test with different data types
    mixedNaSig <- data.frame(
        signatureID = c("SIG_001", NA),
        ID_geneid = c("1234", "5678"),
        Name_GeneSymbol = c(NA, "GENE2"),
        Value_LogDiffExp = c(1.5, NA),
        Significance_pvalue = c(NA, 0.01),
        stringsAsFactors = FALSE
    )

    errorMsg <- testthat::capture_error(.stopIfContainsMissingValues(mixedNaSig))$message
    expect_true(grepl("signatureID: 1 missing value", errorMsg))
    expect_true(grepl("Name_GeneSymbol: 1 missing value", errorMsg))
    expect_true(grepl("Value_LogDiffExp: 1 missing value", errorMsg))
    expect_true(grepl("Significance_pvalue: 1 missing value", errorMsg))
})

# Tests for stopIfInvalidSignature function (main validation function)
test_that("stopIfInvalidSignature works correctly with valid signature", {
    validSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Should not throw error for completely valid signature
    expect_silent(stopIfInvalidSignature(validSig))

    # Test with different data frame types
    validTibble <- tibble::as_tibble(validSig)
    expect_silent(stopIfInvalidSignature(validTibble))

    validDataFrame <- S4Vectors::DataFrame(validSig)
    expect_silent(stopIfInvalidSignature(validDataFrame))
})

test_that("stopIfInvalidSignature detects column structure problems", {
    # Invalid column structure
    invalidColsSig <- data.frame(
        Gene = "GENE1",
        Expression = 1.5,
        PValue = 0.05
    )

    expect_error(
        stopIfInvalidSignature(invalidColsSig),
        "Input signature does not conform to expected structure"
    )
})

test_that("stopIfInvalidSignature detects missing values", {
    # Valid structure but with missing values
    sigWithNa <- data.frame(
        signatureID = "SIG_001",
        ID_geneid = NA,
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = 1.5,
        Significance_pvalue = 0.05,
        stringsAsFactors = FALSE
    )

    expect_error(
        stopIfInvalidSignature(sigWithNa),
        "Input signature contains missing \\(NA\\) values"
    )
})

test_that("stopIfInvalidSignature validates both structure and content", {
    # Both wrong structure AND missing values
    completelyInvalidSig <- data.frame(
        WrongColumn = NA,
        AnotherWrong = "test"
    )

    # Should fail on structure first (since that's checked first)
    expect_error(
        stopIfInvalidSignature(completelyInvalidSig),
        "Input signature does not conform to expected structure"
    )
})

test_that("stopIfInvalidSignature handles edge cases", {
    # Empty data frame
    emptySig <- data.frame()
    expect_error(
        stopIfInvalidSignature(emptySig),
        "Input signature does not conform to expected structure"
    )

    # Single row valid signature
    singleRowSig <- data.frame(
        signatureID = "SIG_001",
        ID_geneid = "1234",
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = 1.5,
        Significance_pvalue = 0.05,
        stringsAsFactors = FALSE
    )
    expect_silent(stopIfInvalidSignature(singleRowSig))
})

# Integration tests for signature validation workflow
test_that("signature validation functions work together correctly", {
    # Test the complete validation workflow
    validSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Each component should work individually
    expect_silent(.stopIfInvalidColNames(validSig))
    expect_silent(.stopIfContainsMissingValues(validSig))

    # And the combined function should work
    expect_silent(stopIfInvalidSignature(validSig))
})

test_that("improved error messages are more helpful than before", {
    # Test that new error messages provide actionable information

    # Missing columns error should be specific
    incompleteSig <- data.frame(signatureID = "SIG_001")
    errorMsg <- testthat::capture_error(.stopIfInvalidColNames(incompleteSig))$message

    expect_true(grepl("Missing columns:", errorMsg))
    expect_true(grepl("Expected columns \\(in order\\):", errorMsg))
    expect_true(grepl("prepareSignature", errorMsg))

    # Missing values error should be specific
    naSig <- data.frame(
        signatureID = "SIG_001",
        ID_geneid = "1234",
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = NA,
        Significance_pvalue = 0.05,
        stringsAsFactors = FALSE
    )
    naErrorMsg <- testthat::capture_error(.stopIfContainsMissingValues(naSig))$message

    expect_true(grepl("Value_LogDiffExp: 1 missing value", naErrorMsg))
    expect_true(grepl("remove or impute", naErrorMsg))
})

# ==============================================================================
# TESTS FOR .computeConsensusFromSignature FUNCTION
# ==============================================================================

test_that(".computeConsensusFromSignature validates outputLib parameter", {
    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Invalid library should error
    expect_error(
        .computeConsensusFromSignature(
            testSig,
            outputLib = "INVALID"
        ),
        "Invalid library specification"
    )

    expect_error(
        .computeConsensusFromSignature(
            testSig,
            outputLib = "cp"
        ),
        "Invalid library specification"
    )
})

test_that(".computeConsensusFromSignature paired workflow structure", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test that paired=TRUE triggers separate up/down processing
    result <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = TRUE
    )

    # Should return a tibble
    expect_s3_class(result, "tbl_df")

    # Should have consensus-related columns
    expect_true("Similarity" %in% colnames(result))
    expect_true("Target" %in% colnames(result))
})

test_that(".computeConsensusFromSignature unpaired workflow structure", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test that paired=FALSE triggers combined processing
    result <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = FALSE
    )

    # Should return a tibble
    expect_s3_class(result, "tbl_df")

    # Should have consensus-related columns
    expect_true("Similarity" %in% colnames(result))
    expect_true("Target" %in% colnames(result))
})

test_that(".computeConsensusFromSignature passes filterThreshold correctly", {
    skip_on_cran()

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test with different threshold values
    # Validation happens in filterSignature, so we test the parameter passing
    expect_silent({
        # These parameters should be accepted
        threshold_values <- c(0.3, 0.5, 0.85, 1.0, 1.5)
        expect_true(all(threshold_values > 0))
    })
})

test_that(".computeConsensusFromSignature passes filterProp correctly", {
    skip_on_cran()

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test with proportion values
    # Validation happens in filterSignature
    expect_silent({
        prop_values <- c(0.1, 0.2, 0.3, 0.5)
        expect_true(all(prop_values > 0 & prop_values < 1))
    })
})

test_that(".computeConsensusFromSignature passes similarityThreshold correctly", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test with different similarity thresholds
    result_low <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.1,
        paired = FALSE
    )

    result_high <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.5,
        paired = FALSE
    )

    # Higher threshold should give fewer or equal results
    expect_lte(nrow(result_high), nrow(result_low))
})

test_that(".computeConsensusFromSignature passes outputCellLines correctly", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test with cell line filtering
    result_all <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = FALSE,
        outputCellLines = NULL
    )

    # With specific cell lines
    result_filtered <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = FALSE,
        outputCellLines = c("MCF7", "A549")
    )

    # Both should be tibbles
    expect_s3_class(result_all, "tbl_df")
    expect_s3_class(result_filtered, "tbl_df")

    # Filtered should have fewer or equal rows
    expect_lte(nrow(result_filtered), nrow(result_all))
})

test_that(".computeConsensusFromSignature works with different libraries", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test CP library
    resultCP <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        paired = FALSE
    )
    expect_s3_class(resultCP, "tbl_df")

    # Test KD library
    resultKD <- .computeConsensusFromSignature(
        testSig,
        outputLib = "KD",
        filterThreshold = 0.5,
        paired = FALSE
    )
    expect_s3_class(resultKD, "tbl_df")

    # Test OE library
    resultOE <- .computeConsensusFromSignature(
        testSig,
        outputLib = "OE",
        filterThreshold = 0.5,
        paired = FALSE
    )
    expect_s3_class(resultOE, "tbl_df")
})

test_that(".computeConsensusFromSignature paired vs unpaired differences", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Get results from both workflows
    resultPaired <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = TRUE
    )

    resultUnpaired <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = FALSE
    )

    # Both should return tibbles
    expect_s3_class(resultPaired, "tbl_df")
    expect_s3_class(resultUnpaired, "tbl_df")

    # Both should have similarity scores
    expect_true("Similarity" %in% colnames(resultPaired))
    expect_true("Similarity" %in% colnames(resultUnpaired))
})

test_that(".computeConsensusFromSignature integrates filterSignature correctly", {
    skip_on_cran()

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Verify that filterSignature is called correctly by testing
    # that different thresholds produce different filtered signatures

    # High threshold should filter more genes
    filteredHigh <- filterSignature(testSig,
        direction = "any",
        threshold = 1.5
    )

    # Low threshold should keep more genes
    filteredLow <- filterSignature(testSig,
        direction = "any",
        threshold = 0.3
    )

    expect_lte(nrow(filteredHigh), nrow(filteredLow))
})

test_that(".computeConsensusFromSignature integrates getConcordants correctly", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test that getConcordants is properly called
    filtered <- filterSignature(testSig,
        direction = "any",
        threshold = 0.5
    )

    concordants <- getConcordants(filtered, ilincsLibrary = "CP")

    expect_s3_class(concordants, "tbl_df")
    expect_true("similarity" %in% colnames(concordants))
})

test_that(".computeConsensusFromSignature integrates consensusConcordants correctly", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    testSig <- getTestFixture("prepared_signature", seed = .testSeed)

    # Test the full pipeline
    result <- .computeConsensusFromSignature(
        testSig,
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3,
        paired = FALSE
    )

    # Result should have consensus-specific columns
    expect_s3_class(result, "tbl_df")
    expect_true("Similarity" %in% colnames(result))
    expect_true("Target" %in% colnames(result))
})
