## Basic seed management -------------------------------------------------------
.testSeed <- 12345L
setTestSeed <- function(seed = .testSeed) {
    set.seed(seed)
    invisible(seed)
}

createInputDge <- function(nGenes = NULL, seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)
    nGenes <- if (!is.null(nGenes)) nGenes else sample(5L:10L, size = 1L)

    # Prefer sampling from package dataset `l1000` to satisfy membership checks
    pool <- tryCatch(
        {
            if (exists("l1000", inherits = TRUE) && !is.null(l1000[["L1000"]])) unique(l1000[["L1000"]]) else NULL
        },
        error = function(e) NULL
    )
    if (is.null(pool) || length(pool) < nGenes) {
        pool <- c("CCND3", "FEZ2", "MAPK9", "TERT", "ANKRD10", "BECN1", "CANT1", "ACAA1", "PXN", "STAT1")
    }
    genes <- sample(pool, nGenes, replace = FALSE)

    tibble::tibble(
        Gene = genes,
        LogFC = stats::rnorm(nGenes, mean = 0L, sd = 1.5),
        PValue = stats::runif(nGenes, min = 1e-4, max = 5e-2)
    )
}

createPreparedDge <- function(nGenes = NULL, seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)
    nGenes <- if (!is.null(nGenes)) nGenes else sample(5L:10L, size = 1L)

    # Prefer sampling from package dataset `l1000` to satisfy membership checks
    pool <- tryCatch(
        {
            if (exists("l1000", inherits = TRUE) && !is.null(l1000[["L1000"]])) unique(l1000[["L1000"]]) else NULL
        },
        error = function(e) NULL
    )
    if (is.null(pool) || length(pool) < nGenes) {
        pool <- c("CCND3", "FEZ2", "MAPK9", "TERT", "ANKRD10", "BECN1", "CANT1", "ACAA1", "PXN", "STAT1")
    }
    genes <- sample(pool, nGenes, replace = FALSE)

    # Map gene symbols -> ENTREZ IDs from l1000 (order-preserving)
    idMap <- tryCatch(
        {
            if (exists("l1000", inherits = TRUE) &&
                !is.null(l1000[["L1000"]]) &&
                !is.null(l1000[["ENTREZID"]])) {
                stats::setNames(as.character(l1000[["ENTREZID"]]), l1000[["L1000"]])
            } else {
                NULL
            }
        },
        error = function(e) NULL
    )

    geneIds <- if (!is.null(idMap)) {
        unname(idMap[genes])
    } else {
        rep(NA_character_, length(genes))
    } # nolint

    tibble::tibble(
        signatureID = rep("LINCSID_00000", length(genes)),
        ID_geneid = geneIds,
        Name_GeneSymbol = genes,
        Value_LogDiffExp = stats::rnorm(nGenes, mean = 0L, sd = 1.5),
        Significance_pvalue = stats::runif(nGenes, min = 1e-4, max = 5e-2)
    )
}

createConcordantsTable <- function(nRows = NULL, library = "CP", seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)
    nRows <- if (!is.null(nRows)) {
        nRows
    } else {
        sample(5L:10L, size = 1L)
    }

    halfNum <- floor(nRows / 2L)

    drugIds <- stringr::str_glue("DRUG_{id}",
        id = stringr::str_pad(1L:nRows, 5L, pad = "0")
    )
    treatmentIds <- stringr::str_glue("GENE_{id}", id = stringr::str_pad(1L:nRows, 3L, pad = "0"))

    compound <- sample(drugIds, size = nRows, replace = TRUE)
    treatment <- sample(treatmentIds, size = nRows, replace = TRUE)
    lincsPertId <- stringr::str_replace(compound, stringr::fixed("DRUG"), "PERT")
    geneTargets <- rep(NA_character_, nRows)
    nGenes <- rep(978L, nRows)
    similarity <- stats::runif(nRows, min = -1L, max = 1L)
    pValue <- stats::runif(nGenes, min = 1e-4, max = 5e-2)
    signatureId <- stringr::str_glue(
        "LINCSID_{id}",
        id = stringr::str_pad(1L:nRows, 5L, pad = "0")
    )
    rowCol <- signatureId
    concentration <- sample(c("10uM", "0.5uM", "0.12uM"), size = nRows, replace = TRUE)
    cellLine <- sample(c("A375", "A549", "HELA", "PC3"), size = nRows, replace = TRUE)
    time <- sample(c("6H", "12H", "24H"), size = nRows, replace = TRUE)


    if (library == "CP") {
        concordanceData <- tibble(
            similarity = similarity, pValue = pValue,
            nGenes = nGenes, compound = compound, lincsPertID = lincsPertId,
            GeneTargets = geneTargets, concentration = concentration,
            time = time, `_row` = rowCol, signatureid = signatureId, cellline = cellLine
        )
    } else {
        concordanceData <- tibble(
            similarity = similarity, pValue = pValue,
            nGenes = nGenes, treatment = treatment, lincsPertID = lincsPertId,
            GeneTargets = geneTargets,
            time = time, `_row` = rowCol, signatureid = signatureId, cellline = cellLine
        )
    }
}

createValidIlincsResponse <- function(nRows = NULL, library = "CP", seed = NULL) {
    concordanceData <- createConcordantsTable(nRows = nRows, library = library, seed = seed) |>
        purrr::pmap(function(...) list(...))

    responseData <- list(status = list(
        sessionID = list(Sys.Date() |> as.character()),
        gpgene = list("NA"),
        gpprobe = list("NA"),
        Remark = list("Done"),
        NoOfGenes = list(978L),
        NoOfProbes = list("NA"),
        concordanceTable = concordanceData,
        corTablePath = list("Empty Path String")
    )) |>
        jsonlite::toJSON(auto_unbox = TRUE)

    httr2::new_response(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5",
        status_code = 200L,
        method = "POST",
        headers = list(`Content-Type` = "application/json"),
        body = charToRaw(responseData)
    )
}

createEmptyIlincsResponse <- function() {
    responseData <- list(status = list(
        sessionID = list(Sys.Date() |> as.character()),
        gpgene = list("NA"),
        gpprobe = list("NA"),
        Remark = list("Done"),
        NoOfGenes = list(978L),
        NoOfProbes = list("NA"),
        concordanceTable = list()
    )) |>
        jsonlite::toJSON(auto_unbox = TRUE)

    httr2::new_response(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5",
        status_code = 200L,
        method = "POST",
        headers = list(`Content-Type` = "application/json"), # nolint: nonportable_path_linter.
        body = charToRaw(responseData)
    )
}

createIlincsErrorResponse400 <- function() {
    responseData <- list(error = list(
        statusCode = 400L,
        name = "An Example Error for testing",
        message = "All your base are belong to us.\n\nYou have no chance to survive make your time.",
        code = "ilincsR_ERROR"
    )) |>
        jsonlite::toJSON(auto_unbox = TRUE)

    httr2::new_response(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5",
        status_code = 400L,
        method = "POST",
        headers = list(`Content-Type` = "application/json"), # nolint: nonportable_path_linter.
        body = charToRaw(responseData)
    )
}

createIlincsErrorResponse500 <- function() {
    responseData <- list(error = list(
        statusCode = 500L,
        message = "Internal Server Error"
    )) |>
        jsonlite::toJSON(auto_unbox = TRUE)

    httr2::new_response(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5",
        status_code = 500L,
        method = "POST",
        headers = list(`Content-Type` = "application/json"), # nolint: nonportable_path_linter.
        body = charToRaw(responseData)
    )
}

## Public test fixture accessor (signature -> DGE-like for these tests) -------
getTestFixture <- function(
  type = c(
      "prepared_signature", "input_signature",
      "valid_ilincs_response", "empty_ilincs_response",
      "error_ilincs_response_400", "error_ilincs_response_500",
      "concordants_table"
  ), seed = NULL, ...
) {
    type <- match.arg(type)
    switch(type,
        input_signature = createInputDge(seed = seed, ...),
        prepared_signature = createPreparedDge(seed = seed, ...),
        valid_ilincs_response = createValidIlincsResponse(seed = seed, ...),
        empty_ilincs_response = createEmptyIlincsResponse(),
        error_ilincs_response_400 = createIlincsErrorResponse400(),
        error_ilincs_response_500 = createIlincsErrorResponse500(),
        concordants_table = createConcordantsTable(seed = seed, ...),
        stop("Unsupported fixture type: ", type, call. = FALSE)
    )
}

################################################################################
# (Optional) Small edge-case helpers retained (lightweight) --------------------
createZeroSignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = "GENE_ZERO",
        Value_LogDiffExp = 0L
    )
}
createSmallPositiveSignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = c("GENE_A", "GENE_B"),
        Value_LogDiffExp = c(0.01, 0.02)
    )
}
createSmallNegativeSignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = c("GENE_A", "GENE_B"),
        Value_LogDiffExp = c(-0.01, -0.02)
    )
}

createEmptySignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = character(),
        Value_LogDiffExp = double(),
        Significance_pvalue = double()
    )
}

lincsKdId <- function() "LINCSKD_13548"
lincsOeId <- function() "LINCSOE_10751"
lincsCpId <- function() "LINCSCP_174580"
