# Test prepareSignature

## Test Invalid Inputs

test_that("prepareSignature throws an error if geneColumn is not present", {
    expect_error(
        prepareSignature(exampleSignature(),
            geneColumn = "Gene",
            logfcColumn = "Value_LogDiffExp",
            pvalColumn = "Significance_pvalue"
        ),
        "geneColumn should be present in the dataframe"
    )
})

test_that("prepareSignature throws an error if logfcColumn is not present", {
    expect_error(
        prepareSignature(exampleSignature(),
            geneColumn = "Name_GeneSymbol",
            logfcColumn = "logFC2",
            pvalColumn = "Significance_pvalue"
        ),
        "logfcColumn should be present in the dataframe"
    )
})

test_that("prepareSignature throws an error if pvalColumn is not present", {
    expect_error(
        prepareSignature(exampleSignature(),
            geneColumn = "Name_GeneSymbol",
            logfcColumn = "Value_LogDiffExp",
            pvalColumn = "PValue2"
        ),
        "pvalColumn should be present in the dataframe"
    )
})

## Test Valid Inputs With PValue

test_that("prepareSignature returns a dataframe with the correct columns", {
    signature <- prepareSignature(exampleSignature(),
        geneColumn = "Name_GeneSymbol",
        logfcColumn = "Value_LogDiffExp",
        pvalColumn = "Significance_pvalue"
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

test_that(
    "prepareSignature returns a dataframe with the correct number of rows",
    {
        signature <- prepareSignature(exampleSignature(),
            geneColumn = "Name_GeneSymbol",
            logfcColumn = "Value_LogDiffExp",
            pvalColumn = "Significance_pvalue"
        )
        expect_lte(nrow(signature), 978L)
    }
)

test_that(
    "prepareSignature returns a dataframe with the correct gene symbols",
    {
        signature <- prepareSignature(exampleSignature(),
            geneColumn = "Name_GeneSymbol",
            logfcColumn = "Value_LogDiffExp",
            pvalColumn = "Significance_pvalue"
        )
        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)


## Test Valid Inputs Without PValue

test_that("prepareSignature returns a dataframe with the correct columns", {
    signature <- prepareSignature(exampleSignature(),
        geneColumn = "Name_GeneSymbol",
        logfcColumn = "Value_LogDiffExp", pvalColumn = NA
    )
    expect_named(
        signature,
        c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
    )
})

test_that(
    "prepareSignature returns a dataframe with the correct number of rows",
    {
        signature <- prepareSignature(exampleSignature(),
            geneColumn = "Name_GeneSymbol",
            logfcColumn = "Value_LogDiffExp",
            pvalColumn = NA
        )
        expect_lte(nrow(signature), 978L)
    }
)

test_that(
    "prepareSignature returns a dataframe with the correct gene symbols",
    {
        signature <- prepareSignature(exampleSignature(),
            geneColumn = "Name_GeneSymbol",
            logfcColumn = "Value_LogDiffExp",
            pvalColumn = NA
        )
        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)
