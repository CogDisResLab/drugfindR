test_that("drugfindRCoreData slot accessors work as expected", {
    # Create a mock instance of drugfindRCoreData
    obj <- new("drugfindRCoreData",
        signature = exampleSignature(),
        filteredSignature = {
            exampleSignature() |> filter(abs(Value_LogDiffExp) > 1)
        },
        filterThresholdUp = 1.0,
        filterThresholdDown = -1.0,
        concordanceLimitUp = 0.5,
        concordanceLimitDown = -0.5,
        unfilteredConcordants = concordantsCp(),
        filteredConcordants = consensusConcordantsCpPaired()
    )

    ## Signature
    expect_equal(getCoreSignature(obj), exampleSignature())

    ## Filtered Signature
    expect_equal(getFilteredSignature(obj), {
        exampleSignature() |> filter(abs(Value_LogDiffExp) > 1)
    })

    ## Filter Thresholds
    expect_equal(getFilterThresholdUp(obj), 1.0)
    expect_equal(getFilterThresholdDown(obj), -1.0)
    expect_equal(getFilterThresholds(obj), c(-1.0, 1.0))

    setFilterThresholdUp(obj) <- 2.0
    expect_equal(getFilterThresholdUp(obj), 2.0)

    setFilterThresholdDown(obj) <- -2.0
    expect_equal(getFilterThresholdDown(obj), -2.0)

    setFilterThresholds(obj) <- 3.0
    expect_equal(getFilterThresholds(obj), c(-3.0, 3.0))

    setFilterThresholds(obj) <- c(-4.0, 4.0)
    expect_equal(getFilterThresholds(obj), c(-4.0, 4.0))

    expect_error(setFilterThresholds(obj) <- c(1, 2, 3), "Incorrect number of items in length")

    ## Concordance Limits
    expect_equal(getConcordanceLimitUp(obj), 0.5)
    expect_equal(getConcordanceLimitDown(obj), -0.5)
    expect_equal(getConcordanceLimits(obj), c(-0.5, 0.5))

    setConcordanceLimitUp(obj) <- 1.0
    expect_equal(getConcordanceLimitUp(obj), 1.0)

    setConcordanceLimitDown(obj) <- -1.0
    expect_equal(getConcordanceLimitDown(obj), -1.0)

    setConcordanceLimits(obj) <- 0.85
    expect_equal(getConcordanceLimits(obj), c(-0.85, 0.85))

    setConcordanceLimits(obj) <- c(-0.5, 0.5)
    expect_equal(getConcordanceLimits(obj), c(-0.5, 0.5))

    expect_error(setConcordanceLimits(obj) <- c(1, 2, 3), "Incorrect number of items in length")

    expect_equal(getUnfilteredConcordants(obj), concordantsCp())

    expect_equal(getFilteredConcordants(obj), consensusConcordantsCpPaired())
})
