test_that("drugfindRCoreData constructor works with minimal input", {
    sig <- tibble::tibble(Gene = letters[1:3], Value_LogDiffExp = c(-1, 0, 1))

    core <- drugfindRCoreData(signature = sig)

    expect_s4_class(core, "drugfindRCoreData")
    expect_equal(core@signature, sig)
    expect_equal(core@filterThresholdDown, 0)
    expect_equal(core@filterThresholdUp, 0)
    expect_null(core@filteredSignature)
    expect_equal(core@concordanceLimitDown, -0.2)
    expect_equal(core@concordanceLimitUp, 0.2)
    expect_null(core@unfilteredConcordants)
    expect_null(core@filteredConcordants)
})

test_that("drugfindRCoreData handles full input correctly", {
    sig <- tibble::tibble(Gene = c("A", "B"), Value_LogDiffExp = c(-1.5, 1.2))
    filtered <- tibble::tibble(Gene = "B", Value_LogDiffExp = 1.2)
    concordants <- tibble::tibble(Drug = "Aspirin", Score = 0.8)

    core <- drugfindRCoreData(
        signature = sig,
        filterThresholdUp = 1,
        filterThresholdDown = -1,
        concordanceLimitUp = 0.75,
        concordanceLimitDown = -0.75
    )

    expect_s4_class(core, "drugfindRCoreData")
    expect_equal(core@filterThresholdDown, -1)
    expect_equal(core@filterThresholdUp, 1)
    expect_equal(core@concordanceLimitDown, -0.75)
    expect_equal(core@concordanceLimitUp, 0.75)
    expect_null(core@filteredSignature)
    expect_null(core@unfilteredConcordants)
    expect_null(core@filteredConcordants)
})
