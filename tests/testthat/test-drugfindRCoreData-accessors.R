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
    setCoreSignature(obj) <- emptySignature()
    expect_equal(getCoreSignature(obj), emptySignature())

    ## Filtered Signature
    expect_equal(getFilteredSignature(obj), {
        exampleSignature() |> filter(abs(Value_LogDiffExp) > 1)
    })
    setFilteredSignature(obj) <- emptySignature()
    expect_equal(getFilteredSignature(obj), emptySignature())

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

    setConcordanceLimitUp(obj) <- 1.5
    expect_equal(getConcordanceLimitUp(obj), 1.5)

    setConcordanceLimitDown(obj) <- -1.5
    expect_equal(getConcordanceLimitDown(obj), -1.5)

    setConcordanceLimits(obj) <- 2.0
    expect_equal(getConcordanceLimits(obj), c(-2.0, 2.0))

    setConcordanceLimits(obj) <- c(-2.5, 2.5)
    expect_equal(getConcordanceLimits(obj), c(-2.5, 2.5))

    expect_error(setConcordanceLimits(obj) <- c(1, 2, 3), "Incorrect number of items in length")

    expect_equal(getUnfilteredConcordants(obj), concordantsCp())
    setUnfilteredConcordants(obj) <- concordantsOe()
    expect_equal(getUnfilteredConcordants(obj), concordantsOe())

    expect_equal(getFilteredConcordants(obj), consensusConcordantsCpPaired())
    setFilteredConcordants(obj) <- concordantsOe()
    expect_equal(getFilteredConcordants(obj), concordantsOe())
})
