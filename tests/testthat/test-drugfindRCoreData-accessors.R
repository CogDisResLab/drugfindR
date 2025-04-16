test_that("getters and setters work correctly for all slots", {
    sig <- tibble::tibble(Gene = c("A", "B"), Value_LogDiffExp = c(-1, 1))
    concordants <- tibble::tibble(Drug = "Aspirin", Score = 0.9)
    obj <- drugfindRCoreData(signature = sig)

    # signature
    expect_equal(getCoreSignature(obj), sig)
    new_sig <- tibble::tibble(Gene = "C", Value_LogDiffExp = 2)
    setCoreSignature(obj) <- new_sig
    expect_equal(getCoreSignature(obj), new_sig)

    # thresholds
    setFilterThresholdUp(obj) <- 0.75
    expect_equal(getFilterThresholdUp(obj), 0.75)

    setFilterThresholdDown(obj) <- -0.5
    expect_equal(getFilterThresholdDown(obj), -0.5)

    setFilterThresholds(obj) <- 0.25
    expect_equal(getFilterThresholdUp(obj), 0.25)
    expect_equal(getFilterThresholdDown(obj), -0.25)
    expect_equal(getFilterThresholds(obj), c(-0.25, 0.25))

    # filteredSignature
    filtered <- tibble::tibble(Gene = "C", Value_LogDiffExp = 2)
    setFilteredSignature(obj) <- filtered
    expect_equal(getFilteredSignature(obj), filtered)

    # concordance limits
    setConcordanceLimitUp(obj) <- 0.8
    expect_equal(getConcordanceLimitUp(obj), 0.8)

    setConcordanceLimitDown(obj) <- -0.8
    expect_equal(getConcordanceLimitDown(obj), -0.8)

    setConcordanceLimits(obj) <- 0.25
    expect_equal(getConcordanceLimitUp(obj), 0.25)
    expect_equal(getConcordanceLimitDown(obj), -0.25)
    expect_equal(getConcordanceLimits(obj), c(-0.25, 0.25))

    # concordants
    setUnfilteredConcordants(obj) <- concordants
    expect_equal(getUnfilteredConcordants(obj), concordants)

    setFilteredConcordants(obj) <- concordants
    expect_equal(getFilteredConcordants(obj), concordants)
})
