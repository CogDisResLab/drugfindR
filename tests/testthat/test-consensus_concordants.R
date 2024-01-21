# Test consensus_concordants

## Test Invalid Inputs

test_that(
    "consensus_concordants throws an error if paired\
    is TRUE and only one dataframe is passed",
    {
        expect_error(
            consensus_concordants(concordants_cp(), paired = TRUE)
        )
    }
)

test_that(
    "consensus_concordants throws an error if paired \
    is FALSE and more than one dataframe is passed",
    {
        expect_error(
            consensus_concordants(
                concordants_cp(), concordants_cp(),
                paired = FALSE
            )
        )
    }
)

## Test Unpaired Analysis

test_that("consensus_concordants properly handles similarity threshold", {
    consensus_concordants_result <-
        consensus_concordants(
            concordants_cp(),
            cutoff = 0.321
        )
    expect_identical(
        nrow(consensus_concordants_result), 4290L
    )
    expect_true(
        all(
            abs(consensus_concordants_result[["Similarity"]]) >= 0.321
        )
    )
    expect_named(
        consensus_concordants_result, consensus_concordants_col_names()
    )
})

test_that("consensus_concordants properly handles single cell line filtering", {
    consensus_concordants_result <-
        consensus_concordants(concordants_cp(), cell_line = "A375")
    expect_identical(nrow(consensus_concordants_result), 858L)
    expect_true(
        all(
            consensus_concordants_result[["TargetCellLine"]] == "A375"
        )
    )
    expect_named(
        consensus_concordants_result, consensus_concordants_col_names()
    )
})

test_that("consensus_concordants properly handles single cell line filtering", {
    consensus_concordants_result <-
        consensus_concordants(
            concordants_cp(),
            cell_line = c("A375", "PC3")
        )
    expect_identical(
        nrow(consensus_concordants_result), 1757L
    )
    expect_true(
        all(
            consensus_concordants_result[["TargetCellLine"]] %in% c("A375", "PC3") # nolint: line_length_linter.
        )
    )
    expect_named(
        consensus_concordants_result, consensus_concordants_col_names()
    )
})

## Test Paired Analysis

test_that("consensus_concordants properly handles paired analysis", {
    concordants <- concordants_cp_paired()
    up_concordants <- concordants[[1L]]
    down_concordants <- concordants[[2L]]
    consensus_concordants_result <-
        consensus_concordants(
            up_concordants, down_concordants,
            paired = TRUE
        )
    expect_identical(nrow(consensus_concordants_result), 1076L)
    expect_true(all(abs(consensus_concordants_result[["Similarity"]]) >= 0.321))
    expect_named(
        consensus_concordants_result, consensus_concordants_col_names()
    )
    expect_false(all(duplicated(consensus_concordants_result[["Target"]])))
})

test_that(
    "consensus_concordants properly handles paired analysis with single cell line filtering", # nolint: line_length_linter.
    {
        concordants <- concordants_cp_paired()
        up_concordants <- concordants[[1L]]
        down_concordants <- concordants[[2L]]
        consensus_concordants_result <- consensus_concordants(
            up_concordants, down_concordants,
            paired = TRUE, cell_line = "A375"
        )
        expect_identical(nrow(consensus_concordants_result), 136L)
        expect_true(
            all(consensus_concordants_result[["TargetCellLine"]] == "A375")
        )
        expect_named(
            consensus_concordants_result, consensus_concordants_col_names()
        )
        expect_false(all(duplicated(consensus_concordants_result[["Target"]])))
    }
)

test_that(
    "consensus_concordants properly handles paired analysis with single cell line filtering", # nolint: line_length_linter.
    {
        concordants <- concordants_cp_paired()
        up_concordants <- concordants[[1L]]
        down_concordants <- concordants[[2L]]
        consensus_concordants_result <- consensus_concordants(
            up_concordants, down_concordants,
            paired = TRUE, cell_line = c("A375", "PC3")
        )
        expect_identical(nrow(consensus_concordants_result), 244L)
        expect_true(
            all(
                consensus_concordants_result[["TargetCellLine"]] %in% c("A375", "PC3") # nolint: line_length_linter.
            )
        )
        expect_named(
            consensus_concordants_result, consensus_concordants_col_names()
        )
        expect_false(
            all(duplicated(consensus_concordants_result[["Target"]]))
        )
    }
)

## Test Consensus OE Concordants

test_that("consensus_concordants properly handles OE concordants", {
    consensus_concordants_result <- consensus_concordants(
        concordants_oe(),
        cutoff = 0.4
    )
    expect_identical(nrow(consensus_concordants_result), 706L)
    expect_true(
        all(
            abs(
                consensus_concordants_result[["Similarity"]]
            ) >= 0.4
        )
    )
    expect_named(
        consensus_concordants_result, consensus_concordants_oe_col_names()
    )
    expect_false(all(duplicated(consensus_concordants_result[["Target"]])))
})
