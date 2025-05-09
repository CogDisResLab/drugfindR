# test_that("summarizeDrugfindRCoreData returns all expected slots", {
#     sig <- tibble::tibble(Gene = "A", Value_LogDiffExp = 1)
#     concordants <- tibble::tibble(Drug = "X", Score = 0.9)
#     obj <- drugfindRCoreData(signature = sig)
#
#     setFilterThresholdUp(obj) <- 1
#     setFilterThresholdDown(obj) <- -1
#     setFilteredSignature(obj) <- sig
#     setConcordanceLimitUp(obj) <- 0.8
#     setConcordanceLimitDown(obj) <- -0.8
#     setUnfilteredConcordants(obj) <- concordants
#     setFilteredConcordants(obj) <- concordants
#
#     summary <- summarizeDrugfindRCoreData(obj)
#
#     expect_named(summary, c(
#         "signature", "filterThresholdUp", "filterThresholdDown",
#         "filteredSignature", "concordanceLimitUp", "concordanceLimitDown",
#         "unfilteredConcordants", "filteredConcordants"
#     ))
#
#     expect_equal(summary$signature, sig)
#     expect_equal(summary$filterThresholdUp, 1)
#     expect_equal(summary$concordanceLimitDown, -0.8)
#     expect_equal(summary$filteredConcordants, concordants)
# })
