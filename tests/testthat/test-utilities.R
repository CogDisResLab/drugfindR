# Test utility functions

library(tibble)
library(S4Vectors)

test_that("targetRename works correctly with treatment column", {
    inputNamesWithTreatment <- c(
        "signature",
        "target", "cellLine", "time", "treatment",
        "similarity", "direction", "pvalue"
    )
    expectedNames <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "Similarity", "SignatureDirection", "pValue"
    )

    result <- targetRename(inputNamesWithTreatment)
    expect_identical(result, expectedNames)
})

test_that("targetRename works correctly without treatment column", {
    inputNamesWithoutTreatment <- c(
        "signature",
        "target",
        "cellLine",
        "time",
        "concentration",
        "similarity",
        "direction",
        "pvalue"
    )
    expectedNames <- c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "Similarity",
        "SignatureDirection", "pValue"
    )

    result <- targetRename(inputNamesWithoutTreatment)
    expect_identical(result, expectedNames)
})

test_that(".ilincsBaseUrl returns correct URL", {
    expectedUrl <- "https://www.ilincs.org/api"
    result <- .ilincsBaseUrl()
    expect_identical(result, expectedUrl)
    expect_type(result, "character")
    expect_length(result, 1L)
})

test_that(".validateLibrary works correctly for valid libraries", {
    expect_true(.validateLibrary("CP"))
    expect_true(.validateLibrary("KD"))
    expect_true(.validateLibrary("OE"))
})

test_that(".validateLibrary works correctly for invalid libraries", {
    expect_false(.validateLibrary("INVALID"))
    expect_false(.validateLibrary("cp")) # case sensitive
    expect_false(.validateLibrary("kd"))
    expect_false(.validateLibrary("oe"))
    expect_false(.validateLibrary(""))
})

test_that("validateLibraries works correctly for valid library vectors", {
    expect_true(validateLibraries(c("CP", "KD")))
    expect_true(validateLibraries(c("OE", "CP", "KD")))
    expect_true(validateLibraries("CP"))
    expect_true(validateLibraries(c("OE"))) # nolint: unnecessary_concatenation_linter.
})

test_that("validateLibraries works correctly for invalid library vectors", {
    expect_false(validateLibraries(c("CP", "INVALID")))
    expect_false(validateLibraries(c("cp", "KD")))
    expect_false(validateLibraries(c("CP", "KD", "INVALID")))
    expect_false(validateLibraries("INVALID"))
})

test_that("stopIfInvalidLibraries stops for invalid libraries", {
    expect_error(
        stopIfInvalidLibraries(c("CP", "INVALID")),
        "Both input and output libraries must be one of 'OE', 'KD', 'CP'"
    )
    expect_error(
        stopIfInvalidLibraries("INVALID"),
        "Both input and output libraries must be one of 'OE', 'KD', 'CP'"
    )
})

test_that("stopIfInvalidLibraries doesn't stop for valid libraries", {
    expect_silent(stopIfInvalidLibraries("CP"))
    expect_silent(stopIfInvalidLibraries(c("CP", "KD")))
    expect_silent(stopIfInvalidLibraries(c("OE", "CP", "KD")))
})

test_that("loadMetadata returns correct metadata for each library", {
    # Test OE metadata
    oeResult <- loadMetadata("OE")
    expect_s3_class(oeResult, "tbl_df")
    expect_identical(oeResult, oeMetadata)

    # Test KD metadata
    kdResult <- loadMetadata("KD")
    expect_s3_class(kdResult, "tbl_df")
    expect_identical(kdResult, kdMetadata)

    # Test CP metadata
    cpResult <- loadMetadata("CP")
    expect_s3_class(cpResult, "tbl_df")
    expect_identical(cpResult, cpMetadata)
})

test_that("loadMetadata throws error for invalid library", {
    expect_error(
        loadMetadata("INVALID"),
        "Invalid library"
    )
    expect_error(
        loadMetadata("cp"),
        "Invalid library"
    )
})

test_that(".returnLibrary returns correct library IDs", {
    expect_identical(.returnLibrary("OE"), "LIB_11", ignore_attr = "names")
    expect_identical(.returnLibrary("KD"), "LIB_6", ignore_attr = "names")
    expect_identical(.returnLibrary("CP"), "LIB_5", ignore_attr = "names")
})

test_that(".returnLibrary throws error for invalid library", {
    expect_error(
        .returnLibrary("INVALID"),
        "Both input and output libraries must be one of 'OE', 'KD', 'CP'"
    )
})

test_that(".returnUserAgent returns valid user agent string", {
    userAgent <- .returnUserAgent()
    expect_type(userAgent, "character")
    expect_length(userAgent, 1L)
    expect_true(grepl("drugfindR", userAgent))
    expect_gt(nchar(userAgent), 0L)
})
