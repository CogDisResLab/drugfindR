test_that("drugfindRCoreData constructor works with minimal input", {
    sig <- tibble::tibble(Gene = letters[1:3], Value_LogDiffExp = c(-1, 0, 1))

    core <- drugfindRCoreData(signature = sig)

    expect_s4_class(core, "drugfindRCoreData")
    expect_equal(core@signature, sig)
    expect_true(is.na(core@filterThresholdUp))
    expect_true(is.na(core@filterThresholdDown))
    expect_null(core@filteredSignature)
    expect_true(is.na(core@concordanceLimitUp))
    expect_true(is.na(core@concordanceLimitDown))
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
        filteredSignature = filtered,
        concordanceLimitUp = 0.75,
        concordanceLimitDown = -0.75,
        unfilteredConcordants = concordants,
        filteredConcordants = concordants
    )

    expect_s4_class(core, "drugfindRCoreData")
    expect_equal(core@filteredSignature, filtered)
    expect_equal(core@unfilteredConcordants, concordants)
    expect_equal(core@concordanceLimitUp, 0.75)
})

test_that("drugfindRCoreData errors on non-tibble input", {
    df <- data.frame(Gene = c("A", "B"), Value_LogDiffExp = c(-1, 1))
    expect_error(
        drugfindRCoreData(signature = df),
        "`signature` must be a tibble"
    )
})
