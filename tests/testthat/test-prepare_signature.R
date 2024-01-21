# Test prepare_signature

## Test Invalid Inputs

test_that("prepare_signature throws an error if gene_column is not present", {
    expect_error(
        prepare_signature(example_signature(),
            gene_column = "Gene",
            logfc_column = "Value_LogDiffExp",
            pval_column = "Significance_pvalue"
        ),
        "gene_column should be present in the dataframe"
    )
})

test_that("prepare_signature throws an error if logfc_column is not present", {
    expect_error(
        prepare_signature(example_signature(),
            gene_column = "Name_GeneSymbol",
            logfc_column = "logFC2",
            pval_column = "Significance_pvalue"
        ),
        "logfc_column should be present in the dataframe"
    )
})

test_that("prepare_signature throws an error if pval_column is not present", {
    expect_error(
        prepare_signature(example_signature(),
            gene_column = "Name_GeneSymbol",
            logfc_column = "Value_LogDiffExp",
            pval_column = "PValue2"
        ),
        "pval_column should be present in the dataframe"
    )
})

## Test Valid Inputs With PValue

test_that("prepare_signature returns a dataframe with the correct columns", {
    signature <- prepare_signature(example_signature(),
        gene_column = "Name_GeneSymbol",
        logfc_column = "Value_LogDiffExp",
        pval_column = "Significance_pvalue"
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
    "prepare_signature returns a dataframe with the correct number of rows",
    {
        signature <- prepare_signature(example_signature(),
            gene_column = "Name_GeneSymbol",
            logfc_column = "Value_LogDiffExp",
            pval_column = "Significance_pvalue"
        )
        expect_lte(nrow(signature), 978L)
    }
)

test_that(
    "prepare_signature returns a dataframe with the correct gene symbols",
    {
        signature <- prepare_signature(example_signature(),
            gene_column = "Name_GeneSymbol",
            logfc_column = "Value_LogDiffExp",
            pval_column = "Significance_pvalue"
        )
        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)


## Test Valid Inputs Without PValue

test_that("prepare_signature returns a dataframe with the correct columns", {
    signature <- prepare_signature(example_signature(),
        gene_column = "Name_GeneSymbol",
        logfc_column = "Value_LogDiffExp", pval_column = NA
    )
    expect_named(
        signature,
        c("signatureID", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp")
    )
})

test_that(
    "prepare_signature returns a dataframe with the correct number of rows",
    {
        signature <- prepare_signature(example_signature(),
            gene_column = "Name_GeneSymbol",
            logfc_column = "Value_LogDiffExp",
            pval_column = NA
        )
        expect_lte(nrow(signature), 978L)
    }
)

test_that(
    "prepare_signature returns a dataframe with the correct gene symbols",
    {
        signature <- prepare_signature(example_signature(),
            gene_column = "Name_GeneSymbol",
            logfc_column = "Value_LogDiffExp",
            pval_column = NA
        )
        expect_true(all(signature[["Name_GeneSymbol"]] %in% l1000[["L1000"]]))
    }
)
