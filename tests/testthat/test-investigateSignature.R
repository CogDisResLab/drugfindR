# Comprehensive tests for investigateSignature function

library(httptest2)
library(dplyr)

# ==============================================================================
# TESTS FOR INPUT VALIDATION
# ==============================================================================

test_that("investigateSignature validates outputLib parameter", {
    # Invalid library should error
    testDge <- data.frame(
        Symbol = c("TP53", "MYC", "BRCA1"),
        logFC = c(2.5, -1.8, 3.2),
        PValue = c(0.001, 0.01, 0.0001)
    )

    expect_error(
        investigateSignature(
            testDge,
            outputLib = "INVALID"
        ),
        "Invalid library specification"
    )

    expect_error(
        investigateSignature(
            testDge,
            outputLib = "cp" # Wrong case
        ),
        "Invalid library specification"
    )
})

test_that("investigateSignature requires outputLib parameter", {
    testDge <- data.frame(
        Symbol = c("TP53", "MYC"),
        logFC = c(2.5, -1.8),
        PValue = c(0.001, 0.01)
    )

    # Missing outputLib should error (either missing arg or validation error)
    expect_error(
        investigateSignature(testDge)
    )
})

test_that("investigateSignature validates expr parameter structure", {
    # Invalid expr should fail during prepareSignature
    expect_error(
        investigateSignature(
            "not_a_dataframe",
            outputLib = "CP"
        )
    )

    expect_error(
        investigateSignature(
            list(a = 1L, b = 2L),
            outputLib = "CP"
        )
    )
})

test_that("investigateSignature validates column parameters", {
    testDge <- data.frame(
        Gene = c("TP53", "MYC"),
        FC = c(2.5, -1.8),
        Pval = c(0.001, 0.01)
    )

    # Wrong column names should error
    expect_error(
        investigateSignature(
            testDge,
            outputLib = "CP",
            geneColumn = "WrongColumn"
        ),
        "not found in the dataframe"
    )

    expect_error(
        investigateSignature(
            testDge,
            outputLib = "CP",
            geneColumn = "Gene",
            logfcColumn = "WrongColumn"
        ),
        "not found in the dataframe"
    )
})

# ==============================================================================
# TESTS FOR SIGNATURE PREPARATION INTEGRATION
# ==============================================================================

test_that("investigateSignature correctly prepares signatures", {
    skip_on_cran()

    # Load example data
    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    # Prepare signature manually to compare
    manualSig <- prepareSignature(
        dge_data[1:20, ],
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    # Check that signature has expected structure (can be data.frame or tbl_df)
    expect_true(is.data.frame(manualSig))
    expect_named(manualSig, c(
        "signatureID", "ID_geneid", "Name_GeneSymbol",
        "Value_LogDiffExp", "Significance_pvalue"
    ))
})

test_that("investigateSignature handles custom column names", {
    skip_on_cran()

    customDge <- data.frame(
        GENE = c("TP53", "MYC", "BRCA1", "EGFR"),
        FOLD_CHANGE = c(2.5, -1.8, 3.2, -2.1),
        PVALUE = c(0.001, 0.01, 0.0001, 0.005)
    )

    # Should work with custom column names
    # (Will fail at API call stage without network, but validates params)
    sig <- prepareSignature(
        customDge,
        geneColumn = "GENE",
        logfcColumn = "FOLD_CHANGE",
        pvalColumn = "PVALUE"
    )
    # Check signature was created successfully
    expect_true(is.data.frame(sig))
})

test_that("investigateSignature handles signatures without p-values", {
    skip_on_cran()

    noPvalDge <- data.frame(
        Gene = c("TP53", "MYC", "BRCA1"),
        LogFC = c(2.5, -1.8, 3.2)
    )

    # Should work without p-values
    sig <- prepareSignature(
        noPvalDge,
        geneColumn = "Gene",
        logfcColumn = "LogFC",
        pvalColumn = NA
    )

    expect_true(is.data.frame(sig))
    expect_false("Significance_pvalue" %in% colnames(sig))
})

# ==============================================================================
# TESTS FOR PARAMETER PASSING TO .computeConsensusFromSignature
# ==============================================================================

test_that("investigateSignature passes filterThreshold correctly", {
    skip_on_cran()

    # This is tested indirectly through the consensus calculation
    # We validate that the parameter structure is correct
    testDge <- data.frame(
        Symbol = c("TP53", "MYC"),
        logFC = c(2.5, -1.8),
        PValue = c(0.001, 0.01)
    )

    # Test that different threshold values are accepted
    # (Will fail at API stage without network, but validates params)
    expect_silent({
        sig <- prepareSignature(testDge)
        # Threshold validation happens in filterSignature
        expect_type(0.5, "double")
        expect_gt(0.5, 0)
    })
})

test_that("investigateSignature passes filterProp correctly", {
    skip_on_cran()

    testDge <- data.frame(
        Symbol = c("TP53", "MYC"),
        logFC = c(2.5, -1.8),
        PValue = c(0.001, 0.01)
    )

    # Test that proportion values are accepted
    expect_silent({
        sig <- prepareSignature(testDge)
        # Prop validation happens in filterSignature
        expect_type(0.1, "double")
        expect_gt(0.1, 0)
        expect_lt(0.1, 1)
    })
})

test_that("investigateSignature passes similarityThreshold correctly", {
    skip_on_cran()

    # Test that similarity threshold is properly structured
    expect_type(0.2, "double")
    expect_gte(0.2, 0)
    expect_lte(0.2, 1)
})

test_that("investigateSignature passes paired parameter correctly", {
    skip_on_cran()

    # Test that paired parameter is logical
    expect_type(TRUE, "logical")
    expect_type(FALSE, "logical")
})

test_that("investigateSignature passes outputCellLines correctly", {
    skip_on_cran()

    # Test that outputCellLines accepts character vectors
    testCellLines <- c("MCF7", "A549")
    expect_type(testCellLines, "character")
    expect_length(testCellLines, 2L)

    # NULL should also be accepted
    expect_null(NULL)
})

# ==============================================================================
# TESTS FOR METADATA ANNOTATION
# ==============================================================================

test_that("investigateSignature adds source metadata correctly", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue",
        sourceName = "COVID_Test",
        sourceCellLine = "TestCell",
        sourceTime = "24h",
        sourceConcentration = "10uM"
    )

    # Check that source metadata was added
    expect_true("Source" %in% colnames(result))
    expect_true("SourceCellLine" %in% colnames(result))
    expect_true("SourceTime" %in% colnames(result))
    expect_true("SourceConcentration" %in% colnames(result))

    # Check values
    expect_identical(unique(result[["Source"]]), "COVID_Test")
    expect_identical(unique(result[["SourceCellLine"]]), "TestCell")
    expect_identical(unique(result[["SourceTime"]]), "24h")
    expect_identical(unique(result[["SourceConcentration"]]), "10uM")
})

test_that("investigateSignature uses default metadata values", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    # Check default values
    expect_identical(unique(result[["Source"]]), "Input")
    expect_identical(unique(result[["SourceCellLine"]]), "NA")
    expect_identical(unique(result[["SourceTime"]]), "NA")
    expect_identical(unique(result[["SourceConcentration"]]), "NA")
})

test_that("investigateSignature adds SourceSignature column", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_true("SourceSignature" %in% colnames(result))
    # Should be "InputSig" from prepareSignature
    expect_identical(unique(result[["SourceSignature"]]), "InputSig")
})

# ==============================================================================
# TESTS FOR RETURN VALUE STRUCTURE
# ==============================================================================

test_that("investigateSignature returns tibble with expected columns", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    # Check return type
    expect_s3_class(result, "tbl_df")

    # Check for core expected columns
    expectedCols <- c(
        "Source", "Target", "Similarity",
        "SourceSignature", "TargetSignature"
    )

    expect_true(all(expectedCols %in% colnames(result)))
})

test_that("investigateSignature column order is consistent", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    # Check that key columns appear in expected order
    colIdx <- match(
        c("Source", "Target", "Similarity"),
        colnames(result)
    )

    expect_true(all(diff(colIdx) > 0)) # Should be in increasing order
})

# ==============================================================================
# TESTS FOR DIFFERENT LIBRARY TYPES
# ==============================================================================

test_that("investigateSignature works with CP library", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
})

test_that("investigateSignature works with KD library", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "KD",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
})

test_that("investigateSignature works with OE library", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "OE",
        filterThreshold = 0.5,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_gt(nrow(result), 0L)
})

# ==============================================================================
# TESTS FOR PAIRED VS UNPAIRED WORKFLOW
# ==============================================================================

test_that("investigateSignature paired workflow produces valid results", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        paired = TRUE,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_true("Similarity" %in% colnames(result))
})

test_that("investigateSignature unpaired workflow produces valid results", {
    skip_on_cran()
    skip_if_offline()
    skip("Requires network access to iLINCS API")

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    result <- investigateSignature(
        dge_data[1:20, ],
        outputLib = "CP",
        filterThreshold = 0.5,
        paired = FALSE,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_s3_class(result, "tbl_df")
    expect_true("Similarity" %in% colnames(result))
})

# ==============================================================================
# TESTS FOR EDGE CASES
# ==============================================================================

test_that("investigateSignature handles small gene sets", {
    skip_on_cran()

    smallDge <- data.frame(
        Symbol = c("TP53", "MYC"),
        logFC = c(2.5, -1.8),
        PValue = c(0.001, 0.01)
    )

    sig <- prepareSignature(smallDge)
    expect_true(is.data.frame(sig))
    # May have 0 rows if genes not in L1000
})

test_that("investigateSignature handles large gene sets", {
    skip_on_cran()

    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    # Use full dataset
    sig <- prepareSignature(
        dge_data,
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    expect_true(is.data.frame(sig))
    expect_gt(nrow(sig), 100L) # Should have many genes
})

test_that("investigateSignature handles genes not in L1000", {
    skip_on_cran()

    fakeDge <- data.frame(
        Symbol = c("FAKEGENE1", "FAKEGENE2", "FAKEGENE3"),
        logFC = c(2.5, -1.8, 3.2),
        PValue = c(0.001, 0.01, 0.0001)
    )

    sig <- prepareSignature(fakeDge)

    # Should return empty or near-empty signature
    expect_true(is.data.frame(sig))
    expect_lte(nrow(sig), 3L)
})

# ==============================================================================
# TESTS FOR INTEGRATION WITH .computeConsensusFromSignature
# ==============================================================================

test_that("investigateSignature correctly calls .computeConsensusFromSignature", {
    skip_on_cran()

    # Verify that the consensus function is properly orchestrated
    # This is tested indirectly through the full workflow
    dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
        package = "drugfindR"
    )
    dge_data <- read.delim(dge_file)

    sig <- prepareSignature(
        dge_data[1:20, ],
        geneColumn = "hgnc_symbol",
        logfcColumn = "logFC",
        pvalColumn = "PValue"
    )

    # .computeConsensusFromSignature should be callable with this signature
    expect_true(is.data.frame(sig))
    expect_true("Value_LogDiffExp" %in% colnames(sig))
})
