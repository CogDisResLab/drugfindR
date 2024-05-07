# Helper functions for testing

# Load packages
library(tibble)

## Create Empty Signature

emptySignature <- function() {
    tibble::tibble(
        signatureID = rep(NA, 978L),
        ID_geneid = rep(NA, 978L),
        Name_GeneSymbol = rep(NA, 978L),
        Value_LogDiffExp = rep(NA, 978L),
        Significance_pvalue = rep(NA, 978L)
    )
}

## Signature Column names

signatureColNames <- function() {
    colnames(emptySignature())
}

## Return an example signature

exampleSignature <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "exampleSignature.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        getSignature("LINCSKD_28") |>
            saveRDS(file = rdsPath)
    }
}

## Generate concordants for a signature

concordantsCp <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "concordantsCp.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        getConcordants(
            {
                exampleSignature() |> filterSignature(threshold = 1.0)
            },
            "CP",
            "any"
        ) |>
            saveRDS(file = rdsPath, compress = "xz")
    }
}

concordantsCpPaired <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "concordantsCpPaired.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        signatureUpregulated <- exampleSignature() |>
            filterSignature(threshold = 1.0, direction = "up")
        signatureDownregulated <- exampleSignature() |>
            filterSignature(threshold = 1.0, direction = "down")
        upConcordants <- getConcordants(
            signatureUpregulated,
            "CP",
            "up"
        )
        downConcordants <- getConcordants(
            signatureDownregulated,
            "CP",
            "down"
        )
        list(upConcordants, downConcordants) |>
            saveRDS(file = rdsPath, compress = "xz")
    }
}

concordantsOe <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "concordantsOe.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        getConcordants(
            {
                exampleSignature() |> filterSignature(threshold = 1.0)
            },
            "OE",
            "any"
        ) |>
            saveRDS(file = rdsPath, compress = "xz")
    }
}

consensusConcordantsCpPaired <- function() {
    rdsPath <- file.path(
        test_path(),
        "fixtures",
        "consensusConcordantsCpPaired.RDS"
    )
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        concordants <- concordantsCpPaired()
        upConcordants <- concordants[[1L]]
        downConcordants <- concordants[[2L]]
        consensusConcordants(upConcordants, downConcordants,
            paired = TRUE
        ) |>
            saveRDS(file = rdsPath, compress = "xz")
    }
}

# Concordants Column Names

concordantsColNames <- function() {
    colnames(concordantsCp())
}

## Consensus CP Concordants Column Names

consensusConcordantsColNames <- function() { # nolint: object_length_linter.
    c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "Similarity", "SignatureDirection"
    )
}

## Consensus OE Concordants Column Names

consensusConcordantsOeColNames <- function() { # nolint: object_length_linter, line_length_linter.
    colnames(consensusConcordants(concordantsOe()))
}
