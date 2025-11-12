# Test prepareSignature and Internal Functions

## Test .validatePrepareSignatureInput

test_that(".validatePrepareSignatureInput validates input types correctly", {
    testData <- getTestFixture("input_signature", seed = .testSeed)

    # Test non-character geneColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, 123L, "Value_LogDiffExp", "Significance_pvalue"
        ),
        "geneColumn must be a single character string"
    )

    # Test non-character logfcColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, "Name_GeneSymbol", 123L, "Significance_pvalue"
        ),
        "logfcColumn must be a single character string"
    )

    # Test invalid pvalColumn type
    expect_error(
        .validatePrepareSignatureInput(
            testData, "Name_GeneSymbol", "Value_LogDiffExp", 123L
        ),
        "pvalColumn must be a character string or NA"
    )

    # Test multiple values for geneColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, c("Gene1", "Gene2"), "Value_LogDiffExp", "Significance_pvalue"
        ),
        "geneColumn must be a single character string"
    )

    # Test multiple values for logfcColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, "Name_GeneSymbol", c("logFC1", "logFC2"), "Significance_pvalue"
        ),
        "logfcColumn must be a single character string"
    )

    # Test multiple values for pvalColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, "Name_GeneSymbol", "Value_LogDiffExp", c("pval1", "pval2")
        ),
        "pvalColumn must be a single character string or NA"
    )
})

test_that(".validatePrepareSignatureInput validates dataframe structure", {
    # Test null dataframe
    expect_error(
        .validatePrepareSignatureInput(
            NULL, "Gene", "LogFC", "PValue"
        ),
        "dge must be a non-empty dataframe-like object"
    )

    # Test empty dataframe
    emptyData <- data.frame()
    expect_error(
        .validatePrepareSignatureInput(
            emptyData, "Gene", "LogFC", "PValue"
        ),
        "dge must be a non-empty dataframe-like object"
    )
})

test_that(".validatePrepareSignatureInput validates column existence", {
    testData <- getTestFixture("input_signature", seed = .testSeed)

    # Test missing geneColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, "NonExistentGene", "LogFC", "PValue"
        ),
        "geneColumn 'NonExistentGene' not found in the dataframe"
    )

    # Test missing logfcColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, "Gene", "NonExistentLogFC", "PValue"
        ),
        "logfcColumn 'NonExistentLogFC' not found in the dataframe"
    )

    # Test missing pvalColumn
    expect_error(
        .validatePrepareSignatureInput(
            testData, "Gene", "LogFC", "NonExistentPvalue"
        ),
        "pvalColumn 'NonExistentPvalue' not found in the dataframe"
    )

    # Test valid call with NA pvalColumn
    expect_silent(
        .validatePrepareSignatureInput(testData, "Gene", "LogFC", NA)
    )

    # Test valid call with all columns present
    expect_silent(
        .validatePrepareSignatureInput(
            testData, "Gene", "LogFC", "PValue"
        )
    )
})

## Test .mapToL1000WithPvalues

test_that(".mapToL1000WithPvalues maps data correctly", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)
    mapped <- .mapToL1000WithPvalues(filteredData, "Gene", "LogFC", "PValue")

    # Check expected columns
    expectedCols <- c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
    expect_named(mapped, expectedCols)

    # Check signatureID is set correctly
    expect_true(all(mapped[["signatureID"]] == "InputSig"))

    # Check gene symbols are from L1000
    expect_true(all(mapped[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))

    # Check data types
    expect_type(mapped[["Value_LogDiffExp"]], "double")
    expect_type(mapped[["Significance_pvalue"]], "double")
})

## Test .mapToL1000WithoutPvalues

test_that(".mapToL1000WithoutPvalues maps data correctly", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed) |>
        select(-PValue)

    mapped <- .mapToL1000WithoutPvalues(filteredData, "Gene", "LogFC")

    # Check expected columns (should not include Significance_pvalue)
    expectedCols <- c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
    expect_named(mapped, expectedCols)

    # Check signatureID is set correctly
    expect_true(all(mapped[["signatureID"]] == "InputSig"))

    # Check gene symbols are from L1000
    expect_true(all(mapped[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))

    # Check data types
    expect_type(mapped[["Value_LogDiffExp"]], "double")
    expect_type(mapped[["Significance_pvalue"]], "NULL")
})

## Test .processToL1000Signature

test_that(".processToL1000Signature dispatches correctly with p-values", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    result <- .processToL1000Signature(filteredData, "Gene", "LogFC", "PValue")

    # Should include p-value column
    expectedCols <- c(
        "signatureID", "ID_geneid",
        "Name_GeneSymbol", "Value_LogDiffExp",
        "Significance_pvalue"
    )
    expect_named(result, expectedCols)
})

test_that(".processToL1000Signature dispatches correctly without p-values", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed) |>
        select(-PValue)

    result <- .processToL1000Signature(filteredData, "Gene", "LogFC", NA)

    # Should not include p-value column in meaningful way
    expectedCols <- c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
    expect_named(result, expectedCols)
})

## Test Main prepareSignature Function - Invalid Inputs

test_that("prepareSignature throws an error if geneColumn is not present", {
    expect_error(
        prepareSignature(getTestFixture("input_signature", seed = .testSeed),
            geneColumn = "Missing",
            logfcColumn = "LogFC",
            pvalColumn = "PValue"
        ),
        "geneColumn 'Missing' not found in the dataframe"
    )
})

test_that("prepareSignature throws an error if logfcColumn is not present", {
    expect_error(
        prepareSignature(getTestFixture("input_signature", seed = .testSeed),
            geneColumn = "Gene",
            logfcColumn = "Missing",
            pvalColumn = "PValue"
        ),
        "logfcColumn 'Missing' not found in the dataframe"
    )
})

test_that("prepareSignature throws an error if pvalColumn is not present", {
    expect_error(
        prepareSignature(getTestFixture("input_signature", seed = .testSeed),
            geneColumn = "Gene",
            logfcColumn = "LogFC",
            pvalColumn = "Missing"
        ),
        "pvalColumn 'Missing' not found in the dataframe"
    )
})

## Test Main prepareSignature Function - Valid Inputs With PValue

test_that("prepareSignature returns a dataframe with the correct columns (with pvalue)", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = "PValue"
    )
    expect_named(
        signature,
        c(
            "signatureID",
            "ID_geneid",
            "Name_GeneSymbol",
            "Value_LogDiffExp",
            "Significance_pvalue"
        )
    )
})

test_that("prepareSignature returns correct data types (with pvalue)", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)
    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = "PValue"
    )

    expect_type(signature[["signatureID"]], "character")
    expect_type(signature[["ID_geneid"]], "character")
    expect_type(signature[["Name_GeneSymbol"]], "character")
    expect_type(signature[["Value_LogDiffExp"]], "double")
    expect_type(signature[["Significance_pvalue"]], "double")
})

test_that("prepareSignature returns correct signatureID (with pvalue)", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)
    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = "PValue"
    )

    expect_true(all(signature[["signatureID"]] == "InputSig"))
})

test_that(
    "prepareSignature returns a dataframe with the correct number of rows (with pvalue)",
    {
        filteredData <- getTestFixture("input_signature", seed = .testSeed)
        signature <- prepareSignature(filteredData,
            geneColumn = "Gene",
            logfcColumn = "LogFC",
            pvalColumn = "PValue"
        )
        expect_identical(nrow(signature), 9L)
    }
)

test_that(
    "prepareSignature returns a dataframe with the correct gene symbols (with pvalue)",
    {
        filteredData <- getTestFixture("input_signature", seed = .testSeed)
        signature <- prepareSignature(filteredData,
            geneColumn = "Gene",
            logfcColumn = "LogFC",
            pvalColumn = "PValue"
        )
        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)

test_that("prepareSignature handles duplicate genes correctly (with pvalue)", {
    skip("Duplicate handling not implemented yet")

    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(testData,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = "PValue"
    )

    # Should handle duplicates (exact behavior depends on implementation)
    expect_gte(nrow(signature), 1L)
    expect_lte(nrow(signature), 3L)
    expect_identical(nrow(signature), 2L)
})

## Test Main prepareSignature Function - Valid Inputs Without PValue

test_that("prepareSignature returns a dataframe with the correct columns (without pvalue)", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC", pvalColumn = NA
    )

    # Check that it has the expected columns but not Significance_pvalue
    expectedCols <- c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
    expect_named(signature, expectedCols)
})

test_that("prepareSignature returns correct data types (without pvalue)", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC", pvalColumn = NA
    )

    expect_type(signature[["signatureID"]], "character")
    expect_type(signature[["ID_geneid"]], "character")
    expect_type(signature[["Name_GeneSymbol"]], "character")
    expect_type(signature[["Value_LogDiffExp"]], "double")
})

test_that("prepareSignature returns correct signatureID (without pvalue)", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC", pvalColumn = NA
    )

    expect_true(all(signature[["signatureID"]] == "InputSig"))
})

test_that(
    "prepareSignature returns a dataframe with the correct number of rows (without pvalue)",
    {
        filteredData <- getTestFixture("input_signature", seed = .testSeed)

        signature <- prepareSignature(filteredData,
            geneColumn = "Gene",
            logfcColumn = "LogFC", pvalColumn = NA
        )

        expect_lte(nrow(signature), 978L)
        expect_gt(nrow(signature), 0L)
        expect_identical(nrow(signature), 9L)
    }
)

test_that(
    "prepareSignature returns a dataframe with the correct gene symbols (without pvalue)",
    {
        filteredData <- getTestFixture("input_signature", seed = .testSeed)

        signature <- prepareSignature(filteredData,
            geneColumn = "Gene",
            logfcColumn = "LogFC", pvalColumn = NA
        )

        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)

test_that("prepareSignature handles duplicate genes correctly (without pvalue)", {
    skip("Duplicate handling not implemented yet")
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC", pvalColumn = NA
    )

    # Should handle duplicates (exact behavior depends on implementation)
    expect_gte(nrow(signature), 1L)
})

## Test Edge Cases

test_that("prepareSignature handles empty filtered results", {
    # Create test data with no L1000 genes
    testData <- data.frame(
        Gene = c("FAKEGENE1", "FAKEGENE2", "FAKEGENE3"),
        LogFC = c(1.5, -2.1, 0.5),
        PValue = c(0.01, 0.05, 0.1)
    )

    signature <- prepareSignature(testData,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = "PValue"
    )

    # Should return empty dataframe with correct structure
    expect_identical(nrow(signature), 0L)
    expectedCols <- c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
    expect_named(signature, expectedCols)
})

test_that("prepareSignature handles single gene input", {
    testData <- data.frame(
        Gene = "TP53",
        LogFC = 1.5,
        PValue = 0.01
    )

    signature <- prepareSignature(testData,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = "PValue"
    )

    # Should successfully process single gene
    expect_identical(nrow(signature), 1L)
    expect_identical(
        signature[["Name_GeneSymbol"]][1L],
        "TP53"
    )
})

test_that("prepareSignature works with different column names", {
    testData <- data.frame(
        GENE_SYMBOL = c("TP53", "MYC", "GAPDH"),
        FOLD_CHANGE = c(1.5, -2.1, 0.8),
        P_VALUE = c(0.01, 0.05, 0.02)
    )

    signature <- prepareSignature(testData,
        geneColumn = "GENE_SYMBOL",
        logfcColumn = "FOLD_CHANGE",
        pvalColumn = "P_VALUE"
    )

    # Should work with custom column names
    expect_gt(nrow(signature), 0L)
    expectedCols <- c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")
    expect_true(all(expectedCols %in% names(signature)))
})


## Test Valid Inputs Without PValue

test_that("prepareSignature returns a dataframe with the correct columns", {
    filteredData <- getTestFixture("input_signature", seed = .testSeed)

    signature <- prepareSignature(filteredData,
        geneColumn = "Gene",
        logfcColumn = "LogFC", pvalColumn = NA
    )

    expect_named(
        signature,
        c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
    )
})

test_that(
    "prepareSignature returns a dataframe with the correct number of rows",
    {
        filteredData <- getTestFixture("input_signature", seed = .testSeed)

        signature <- prepareSignature(filteredData,
            geneColumn = "Gene",
            logfcColumn = "LogFC", pvalColumn = NA
        )

        expect_gte(nrow(signature), 1L)
        expect_lte(nrow(signature), 10L)
        expect_identical(nrow(signature), 9L)
    }
)

test_that(
    "prepareSignature returns a dataframe with the correct gene symbols",
    {
        filteredData <- getTestFixture("input_signature", seed = .testSeed)

        signature <- prepareSignature(filteredData,
            geneColumn = "Gene",
            logfcColumn = "LogFC", pvalColumn = NA
        )

        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)
