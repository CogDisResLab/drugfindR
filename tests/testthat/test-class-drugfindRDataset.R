test_that("Constructor works with valid inputs", {
    sig <- exampleSignature()
    dfc <- new("drugfindRCoreData",
        signature = sig,
        inputLibrary = "OE",
        outputLibrary = "CP",
        pairedAnalysis = TRUE,
        cellLines = c("MCF7", "PC3"),
        filterThresholdUp = 1.0,
        filterThresholdDown = -1.0,
        concordanceLimitUp = 0.5,
        concordanceLimitDown = -0.5
    )

    dfds <- new("drugfindRDataset",
        core = dfc,
        source_class = "DFrame"
    )
    expect_s4_class(dfds, "drugfindRDataset")
    expect_s4_class(dfds@core, "drugfindRCoreData")
    expect_equal(dfds@source_class, "DFrame")
    expect_null(getUnfilteredConcordants(dfds@core))
    expect_null(getFilteredConcordants(dfds@core))
})

test_that("Default prototype values are set correctly", {
    dfds <- new("drugfindRDataset")
    expect_s4_class(dfds, "drugfindRDataset")
    expect_s4_class(dfds@core, "drugfindRCoreData")
    expect_equal(dfds@source_class, "tbl_df")
    expect_null(getUnfilteredConcordants(dfds@core))
    expect_null(getFilteredConcordants(dfds@core))
})
