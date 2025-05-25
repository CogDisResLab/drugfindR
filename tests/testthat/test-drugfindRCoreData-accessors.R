test_that("drugfindRCoreData slot accessors work as expected", {
    # Create a mock instance of drugfindRCoreData
    obj <- new("drugfindRCoreData",
        coreSignature = exampleSignature(),
        filteredSignature = {
            exampleSignature() |> filter(abs(Value_LogDiffExp) > 1L)
        },
        filterThresholdUp = 1.0,
        filterThresholdDown = -1.0,
        concordanceLimitUp = 0.5,
        concordanceLimitDown = -0.5,
        unfilteredConcordants = concordantsCp(),
        filteredConcordants = consensusConcordantsCpPaired()
    )

    ## Signature
    expect_identical(coreSignature(obj), exampleSignature())

    ## Filtered Signature
    expect_identical(filteredSignature(obj), {
        exampleSignature() |> filter(abs(Value_LogDiffExp) > 1L)
    })

    ## Filter Thresholds
    expect_equal(filterThresholdUp(obj), 1.0)
    expect_equal(filterThresholdDown(obj), -1.0)
    expect_equal(filterThresholds(obj), c(-1.0, 1.0))

    filterThresholdUp(obj) <- 2.0
    expect_equal(filterThresholdUp(obj), 2.0)

    filterThresholdDown(obj) <- -2.0
    expect_equal(filterThresholdDown(obj), -2.0)

    filterThresholds(obj) <- 3.0
    expect_equal(filterThresholds(obj), c(-3.0, 3.0))

    filterThresholds(obj) <- c(-4.0, 4.0)
    expect_equal(filterThresholds(obj), c(-4.0, 4.0))

    expect_error(
        filterThresholds(obj) <- c(1L, 2L, 3L),
        "Incorrect number of items in length"
    )

    ## Concordance Limits
    expect_equal(concordanceLimitUp(obj), 0.5)
    expect_equal(concordanceLimitDown(obj), -0.5)
    expect_equal(concordanceLimits(obj), c(-0.5, 0.5))

    concordanceLimitUp(obj) <- 1.0
    expect_equal(concordanceLimitUp(obj), 1.0)

    concordanceLimitDown(obj) <- -1.0
    expect_equal(concordanceLimitDown(obj), -1.0)

    concordanceLimits(obj) <- 0.85
    expect_equal(concordanceLimits(obj), c(-0.85, 0.85))

    concordanceLimits(obj) <- c(-0.5, 0.5)
    expect_equal(concordanceLimits(obj), c(-0.5, 0.5))

    expect_error(
        concordanceLimits(obj) <- c(1L, 2L, 3L),
        "Incorrect number of items in length"
    )

    expect_identical(unfilteredConcordants(obj), concordantsCp())

    expect_identical(
        filteredConcordants(obj),
        consensusConcordantsCpPaired()
    )
})
