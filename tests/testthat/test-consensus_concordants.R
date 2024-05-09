# Test consensusConcordants

## Test Invalid Inputs

test_that(
    "consensusConcordants throws an error if paired\
    is TRUE and only one dataframe is passed",
    {
        expect_error(
            consensusConcordants(concordantsCp(), paired = TRUE)
        )
    }
)

test_that(
    "consensusConcordants throws an error if paired \
    is FALSE and more than one dataframe is passed",
    {
        expect_error(
            consensusConcordants(
                concordantsCp(), concordantsCp(),
                paired = FALSE
            )
        )
    }
)

## Test Unpaired Analysis

test_that("consensusConcordants properly handles similarity threshold", {
    consensusConcordantsResult <-
        consensusConcordants(
            concordantsCp(),
            cutoff = 0.321
        )
    expect_identical(
        nrow(consensusConcordantsResult), 4290L
    )
    expect_true(
        all(
            abs(consensusConcordantsResult[["Similarity"]]) >= 0.321
        )
    )
    expect_named(
        consensusConcordantsResult, consensusConcordantsColNames()
    )
})

test_that("consensusConcordants properly handles single cell line filtering", {
    consensusConcordantsResult <-
        consensusConcordants(concordantsCp(), cellLine = "A375")
    expect_identical(nrow(consensusConcordantsResult), 858L)
    expect_true(
        all(
            consensusConcordantsResult[["TargetCellLine"]] == "A375"
        )
    )
    expect_named(
        consensusConcordantsResult, consensusConcordantsColNames()
    )
})

test_that("consensusConcordants properly handles single cell line filtering", {
    consensusConcordantsResult <-
        consensusConcordants(
            concordantsCp(),
            cellLine = c("A375", "PC3")
        )
    expect_identical(
        nrow(consensusConcordantsResult), 1757L
    )
    expect_true(
        all(
            consensusConcordantsResult[["TargetCellLine"]] %in% c("A375", "PC3") # nolint: line_length_linter.
        )
    )
    expect_named(
        consensusConcordantsResult, consensusConcordantsColNames()
    )
})

## Test Paired Analysis

test_that("consensusConcordants properly handles paired analysis", {
    concordants <- concordantsCpPaired()
    upConcordants <- concordants[[1L]]
    downConcordants <- concordants[[2L]]
    consensusConcordantsResult <-
        consensusConcordants(
            upConcordants, downConcordants,
            paired = TRUE
        )
    expect_identical(nrow(consensusConcordantsResult), 1076L)
    expect_true(all(abs(consensusConcordantsResult[["Similarity"]]) >= 0.321))
    expect_named(
        consensusConcordantsResult, consensusConcordantsColNames()
    )
    expect_false(all(duplicated(consensusConcordantsResult[["Target"]])))
})

test_that(
    "consensusConcordants properly handles paired analysis with single cell line filtering", # nolint: line_length_linter.
    {
        concordants <- concordantsCpPaired()
        upConcordants <- concordants[[1L]]
        downConcordants <- concordants[[2L]]
        consensusConcordantsResult <- consensusConcordants(
            upConcordants, downConcordants,
            paired = TRUE, cellLine = "A375"
        )
        expect_identical(nrow(consensusConcordantsResult), 136L)
        expect_true(
            all(consensusConcordantsResult[["TarCellLine"]] == "A375")
        )
        expect_named(
            consensusConcordantsResult, consensusConcordantsColNames()
        )
        expect_false(all(duplicated(consensusConcordantsResult[["Target"]])))
    }
)

test_that(
    "consensusConcordants properly handles paired analysis with single cell line filtering", # nolint: line_length_linter.
    {
        concordants <- concordantsCpPaired()
        upConcordants <- concordants[[1L]]
        downConcordants <- concordants[[2L]]
        consensusConcordantsResult <- consensusConcordants(
            upConcordants, downConcordants,
            paired = TRUE, cellLine = c("A375", "PC3")
        )
        expect_identical(nrow(consensusConcordantsResult), 244L)
        expect_true(
            all(
                consensusConcordantsResult[["TargetCellLine"]] %in% c("A375", "PC3") # nolint: line_length_linter.
            )
        )
        expect_named(
            consensusConcordantsResult, consensusConcordantsColNames()
        )
        expect_false(
            all(duplicated(consensusConcordantsResult[["Target"]]))
        )
    }
)

## Test Consensus OE Concordants

test_that("consensusConcordants properly handles OE concordants", {
    consensusConcordantsResult <- consensusConcordants(
        concordantsOe(),
        cutoff = 0.4
    )
    expect_identical(nrow(consensusConcordantsResult), 706L)
    expect_true(
        all(
            abs(
                consensusConcordantsResult[["Similarity"]]
            ) >= 0.4
        )
    )
    expect_named(
        consensusConcordantsResult, consensusConcordantsOeColNames()
    )
    expect_false(all(duplicated(consensusConcordantsResult[["Target"]])))
})
