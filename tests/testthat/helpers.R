# Helper functions for testing

# Load packages
library(tibble)
library(httr2)

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

## Mocking responses for httr2
mockResponse <- function(response, status) {
    function(req) {
        # Extract serializable request information
        requestInfo <- list(
            method = req$method %||% "GET",
            url = req$url %||% "unknown",
            headers = req$headers %||% list(),
            options = req$options %||% list()
        )

        httr2::response_json(
            url = "https://www.ilincs.org/api",
            body = list(request = requestInfo, response = response),
            status_code = status
        )
    }
}

## Example Response Data

exampleResponse <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "exampleResponse.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        resp <- exampleSignature() |>
            .prepareSignatureFile() |>
            .generateIlincsRequest(
                "CP"
            ) |>
            .executeIlincsRequest()
        respJson <- resp_body_json(resp)

        respJson[["status"]][["concordanceTable"]] <- head(
            respJson[["status"]][["concordanceTable"]]
        )

        updatedResponse <- response_json(
            url = resp_url(resp),
            method = req_get_method(resp),
            status_code = resp_status(resp),
            body = respJson,
            headers = resp_headers(resp)
        )
        saveRDS(updatedResponse, file = rdsPath)
        updatedResponse
    }
}

## Empty concordance Table Response
emptyResponse <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "emptyResponse.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        zeroSignature <- exampleSignature() |>
            mutate(
                Value_LogDiffExp = 0L,
                Significance_pvalue = 0L
            )

        resp <- zeroSignature |>
            .prepareSignatureFile() |>
            .generateIlincsRequest(
                "CP"
            ) |>
            .executeIlincsRequest()
        saveRDS(resp, file = rdsPath)
        resp
    }
}

## Error Response (400)
errorResponse400 <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "errorResponse400.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        resp <- .generateIlincsRequest(here::here("DESCRIPTION"), "CP") |>
            .executeIlincsRequest()
        saveRDS(resp, file = rdsPath)
        resp
    }
}

## Error Response (500)
errorResponse500 <- function() {
    rdsPath <- file.path(test_path(), "fixtures", "errorResponse500.RDS")
    if (file.exists(rdsPath)) {
        readr::read_rds(rdsPath)
    } else {
        resp <- emptySignature() |>
            .prepareSignatureFile() |>
            .generateIlincsRequest("CP") |>
            .executeIlincsRequest()
        saveRDS(resp, file = rdsPath)
        resp
    }
}

## Example request for testing
exampleRequest <- function() {
    exampleSignature() |>
        .prepareSignatureFile() |>
        .generateIlincsRequest(
            "CP"
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

exampleSignatureUpFilter <- function() {
    exampleSignature() |>
        filter(Value_LogDiffExp >= 1.5)
}

exampleSignatureDownFilter <- function() {
    exampleSignature() |>
        filter(Value_LogDiffExp <= -1.5)
}

exampleSignatureAnyFilter <- function() {
    exampleSignature() |>
        filter(abs(Value_LogDiffExp) >= 1.5)
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
            "CP"
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
            "CP"
        )
        downConcordants <- getConcordants(
            signatureDownregulated,
            "CP"
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
            "OE"
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
        "TargetTime", "TargetConcentration", "Similarity",
        "SignatureDirection", "pValue"
    )
}

## Consensus OE Concordants Column Names

consensusConcordantsOeColNames <- function() { # nolint: object_length_linter, line_length_linter.
    colnames(consensusConcordants(concordantsOe()))
}
