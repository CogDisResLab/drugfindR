# Helper functions for testing
# Restructured to minimize duplication and use dummy data

# Load packages
library(tibble)
library(httr2)
library(dplyr)

# ==============================================================================
# 0. REPRODUCIBILITY SETUP
# ==============================================================================

# Set global test seed for reproducibility
.testSeed <- 12345L

# Helper function to set seed before random operations
setTestSeed <- function(seed = .testSeed) {
    set.seed(seed)
    invisible(seed)
}

# Create reproducible test environment
withTestSeed <- function(expr, seed = .testSeed) {
    oldSeed <- .Random.seed
    on.exit(
        {
            if (exists("oldSeed")) {
                .Random.seed <<- oldSeed # nolint: undesirable_operator_linter, object_name_linter.
            }
        },
        add = TRUE
    )

    setTestSeed(seed)
    expr
}

# Generate deterministic test data for consistent results
createReproducibleData <- function(type = c("signature", "concordants"), ..., seed = .testSeed) {
    type <- match.arg(type)
    setTestSeed(seed)

    switch(type,
        signature = createDummySignature(...),
        concordants = createDummyConcordants(...)
    )
}

# ==============================================================================
# 1. DATA SCHEMA DEFINITIONS
# ==============================================================================

# Column name definitions - single source of truth
signatureColumns <- function() {
    c(
        "signatureID", "ID_geneid", "Name_GeneSymbol",
        "Value_LogDiffExp", "Significance_pvalue", "is_L1000"
    )
}

concordantsColumns <- function(library = "CP") {
    baseCols <- c(
        "similarity", "pValue", "nGenes", "compound",
        "lincsPertID", "GeneTargets", "concentration",
        "time", "_row", "signatureid", "cellline"
    )

    # Add library-specific columns
    switch(library,
        CP = c(baseCols, "sig_type", "sig_direction"),
        KD = c(baseCols, "sig_type", "sig_direction"),
        OE = c(baseCols, "sig_type", "sig_direction"),
        baseCols
    )
}

consensusConcordantsColumns <- function(library = "CP") {
    c(
        "TargetSignature", "Target", "TargetCellLine",
        "TargetTime", "TargetConcentration", "Similarity",
        "SignatureDirection", "pValue"
    )
}

# ==============================================================================
# 2. CORE DATA GENERATION (DUMMY DATA)
# ==============================================================================

# Generate realistic test signature data without API calls
# @param seed Optional seed for reproducible random generation
createDummySignature <- function(nGenes = 978L, sigId = lincsKdId(), isL1000 = TRUE, seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)

    tibble::tibble(
        signatureID = rep(sigId, nGenes),
        ID_geneid = seq_len(nGenes),
        Name_GeneSymbol = paste0("GENE", seq_len(nGenes)),
        Value_LogDiffExp = stats::rnorm(nGenes, mean = 0L, sd = 2L),
        Significance_pvalue = stats::runif(nGenes, min = 0.001, max = 0.05),
        is_L1000 = rep(isL1000, nGenes)
    )
}

# Generate empty signature structure
createEmptySignature <- function(nGenes = 978L) {
    tibble::tibble(
        signatureID = rep(NA_character_, nGenes),
        ID_geneid = rep(NA_integer_, nGenes),
        Name_GeneSymbol = rep(NA_character_, nGenes),
        Value_LogDiffExp = rep(NA_real_, nGenes),
        Significance_pvalue = rep(NA_real_, nGenes),
        is_L1000 = rep(NA, nGenes)
    )
}

# Generate realistic test concordants data
# @param seed Optional seed for reproducible random generation
createDummyConcordants <- function(nEntries = 6L, library = "CP", direction = "Up", seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)

    sigType <- switch(library,
        CP = "Chemical Perturbagen",
        KD = "Gene Knockdown",
        OE = "Gene Overexpression",
        "Unknown"
    )

    baseData <- tibble::tibble(
        similarity = stats::runif(nEntries, min = 0.3, max = 0.8),
        pValue = 10L^stats::runif(nEntries, min = -80L, max = -10L),
        nGenes = rep(978L, nEntries),
        compound = paste0("COMPOUND_", seq_len(nEntries)),
        lincsPertID = paste0("LSM-", seq_len(nEntries)),
        GeneTargets = "NA",
        concentration = rep("10uM", nEntries),
        time = rep("6h", nEntries),
        `_row` = paste0("LINCS", toupper(library), "_", seq_len(nEntries)),
        signatureid = paste0("LINCS", toupper(library), "_", seq_len(nEntries)),
        cellline = rep("PC3", nEntries),
        sig_type = rep(sigType, nEntries),
        sig_direction = rep(direction, nEntries)
    )

    baseData
}

# Generate empty concordants structure
createEmptyConcordants <- function(library = "CP") {
    cols <- concordantsColumns(library)
    emptyData <- setNames(
        vector("list", length(cols)),
        cols
    )
    tibble::as_tibble(emptyData)[0L, ]
}

# ==============================================================================
# 3. UNIFIED HTTP MOCK SYSTEM
# ==============================================================================

# Single configurable mock response generator
mockHttpResponse <- function(endpoint, status = 200L, body = NULL, headers = list()) {
    function(req) {
        baseUrl <- switch(endpoint,
            signature = "https://www.ilincs.org/api/ilincsR/downloadSignature",
            concordants = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze",
            "https://www.ilincs.org/api"
        )

        httr2::response_json(
            url = baseUrl,
            status_code = status,
            body = body %||% list(),
            headers = headers
        )
    }
}

# Mock signature API responses
mockSignatureResponse <- function(sigId = lincsKdId(), status = 200L, nGenes = 2L, seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)

    if (status != 200L) {
        return(mockHttpResponse(
            "signature",
            status = status,
            body = list(error = "Signature not found")
        ))
    }

    # Generate realistic signature response based on actual JSON structure
    signatureData <- lapply(seq_len(nGenes), function(i) {
        list(
            signatureID = sigId,
            PROBE = "NA",
            ID_geneid = i,
            Name_GeneSymbol = paste0("GENE", i),
            Value_LogDiffExp = stats::rnorm(1L, sd = 2L),
            Significance_pvalue = stats::runif(1L, 0.001, 0.05)
        )
    })

    body <- list(
        data = list(
            fileName = paste0("sig_", format(Sys.time(), "%a_%b__%d_%H_%M_%S_%Y"), "_test"),
            signature = signatureData
        )
    )

    mockHttpResponse("signature", status = status, body = body)
}

# Mock concordants API responses
mockConcordantsResponse <- function(library = "CP", status = 200L, nEntries = 6L, seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)

    if (status != 200L) {
        return(mockHttpResponse(
            "concordants",
            status = status,
            body = list(error = "Processing failed")
        ))
    }

    # Generate realistic concordants response based on actual JSON structure
    concordantsData <- lapply(seq_len(nEntries), function(i) {
        list(
            similarity = stats::runif(1L, 0.3, 0.8),
            pValue = 10L^stats::runif(1L, -80L, -10L),
            nGenes = 978L,
            compound = paste0("COMPOUND_", i),
            lincsPertID = paste0("LSM-", i),
            GeneTargets = "NA",
            concentration = "10uM",
            time = "6h",
            `_row` = paste0("LINCS", toupper(library), "_", i),
            signatureid = paste0("LINCS", toupper(library), "_", i),
            cellline = "PC3"
        )
    })

    body <- list(
        status = list(
            sessionID = paste0("test_session_", Sys.time()),
            gpgene = "NA",
            gpprobe = "NA",
            Remark = "Done",
            NoOfGenes = 978L,
            NoOfProbes = "NA",
            concordanceTable = concordantsData,
            corTablePath = "/tmp/test_correlation_table.xls"
        )
    )

    mockHttpResponse("concordants", status = status, body = body)
}

# Generic error response generator
mockErrorResponse <- function(status = 404L, message = "Not found") {
    mockHttpResponse("error",
        status = status,
        body = list(error = message)
    )
}

# Simple mock response for general testing
mockResponse <- function(body = list(), status = 200L) {
    function(req) {
        httr2::response_json(
            url = req$url %||% "https://test.example.com",
            status_code = status,
            body = body
        )
    }
}

# ==============================================================================
# 4. DATA TRANSFORMATION UTILITIES
# ==============================================================================

# Filter signature by direction with configurable threshold
filterSignatureByDirection <- function(signature, direction = c("up", "down", "any"), threshold = 1.5) {
    direction <- match.arg(direction)

    switch(direction,
        up = dplyr::filter(signature, Value_LogDiffExp >= threshold),
        down = dplyr::filter(signature, Value_LogDiffExp <= -threshold),
        any = dplyr::filter(signature, abs(Value_LogDiffExp) >= threshold)
    )
}

# Add metadata to signature
addSignatureMetadata <- function(signature, isL1000 = TRUE, sigId = NULL) {
    if (!is.null(sigId)) {
        signature[["signatureID"]] <- sigId
    }
    signature[["is_L1000"]] <- isL1000
    signature
}

# ==============================================================================
# 5. TEST DATA FACTORIES
# ==============================================================================

# Factory for different signature types
signatureFactory <- function(type = c("L1000_KD", "L1000_OE", "L1000_CP", "custom"), nGenes = 978L, seed = NULL) {
    type <- match.arg(type)

    switch(type,
        L1000_KD = createDummySignature(nGenes, lincsKdId(), TRUE, seed),
        L1000_OE = createDummySignature(nGenes, lincsOeId(), TRUE, seed),
        L1000_CP = createDummySignature(nGenes, lincsCpId(), TRUE, seed),
        custom = createDummySignature(nGenes, "CUSTOM_SIG", FALSE, seed)
    )
}

# Factory for concordants data
concordantsFactory <- function(
    library = c("CP", "KD", "OE"),
    direction = c("Up", "Down", "Any"),
    nEntries = 6L, seed = NULL) {
    library <- match.arg(library)
    direction <- match.arg(direction)
    createDummyConcordants(nEntries, library, direction, seed)
}

# Factory for consensus data
consensusFactory <- function(paired = FALSE, library = "CP", seed = NULL) {
    if (paired) {
        upConcordants <- concordantsFactory(library, "Up", seed = seed)
        downConcordants <- concordantsFactory(library, "Down", seed = seed)
        list(upConcordants, downConcordants)
    } else {
        concordantsFactory(library, "Any", seed = seed)
    }
}

# ==============================================================================
# 6. FIXTURE MANAGEMENT SYSTEM
# ==============================================================================

# Centralized fixture management using dummy data
getTestFixture <- function(
    type = c("signature", "concordants", "consensus"),
    library = NULL, cached = TRUE, seed = NULL) {
    type <- match.arg(type)

    fixtureName <- switch(type,
        signature = "exampleSignature.RDS",
        concordants = paste0("concordants", library %||% "Cp", ".RDS"),
        consensus = paste0("consensus", library %||% "Cp", ".RDS")
    )

    fixturePath <- file.path(test_path(), "fixtures", fixtureName)

    if (cached && file.exists(fixturePath)) {
        readr::read_rds(fixturePath)
    } else {
        # Generate dummy data instead of real API calls
        data <- switch(type,
            signature = signatureFactory("L1000_KD", seed = seed),
            concordants = concordantsFactory(library %||% "CP", seed = seed),
            consensus = consensusFactory(library = library %||% "CP", seed = seed)
        )

        if (cached) {
            dir.create(dirname(fixturePath), recursive = TRUE, showWarnings = FALSE)
            saveRDS(data, fixturePath, compress = "xz")
        }
        data
    }
}

# ==============================================================================
# 7. TESTING UTILITIES AND EXAMPLES
# ==============================================================================

# Example request object for testing
exampleRequest <- function() {
    httr2::request("https://www.ilincs.org") |>
        httr2::req_url_path_append("api") |>
        httr2::req_url_path_append("SignatureMeta") |>
        httr2::req_url_path_append("uploadAndAnalyze")
}

# Example successful response for testing
exampleResponse <- function() {
    # Create a realistic concordants response structure
    concordantsData <- list(
        list(
            signatureid = "CPC_001",
            compound = "CHEMBL1271119",
            concentration = "10uM",
            time = "24h",
            cellline = "PC3",
            similarity = 0.54640,
            pValue = 0.001234567890123456789
        ),
        list(
            signatureid = "CPC_002",
            compound = "MLS002264408",
            concentration = "5uM",
            time = "48h",
            cellline = "MCF7",
            similarity = -0.32180,
            pValue = 0.045678901234567890123
        ),
        list(
            signatureid = "CPC_003",
            compound = "VU0413807-2",
            concentration = "1uM",
            time = "6h",
            cellline = "HeLa",
            similarity = 0.78945,
            pValue = 0.000123456789012345678
        ),
        list(
            signatureid = "CPC_004",
            compound = "Flunisolide",
            concentration = "50uM",
            time = "12h",
            cellline = "A549",
            similarity = -0.65432,
            pValue = 0.023456789012345678901
        ),
        list(
            signatureid = "CPC_005",
            compound = "MLS002264499",
            concentration = "25uM",
            time = "72h",
            cellline = "U2OS",
            similarity = 0.41237,
            pValue = 0.012345678901234567890
        ),
        list(
            signatureid = "CPC_006",
            compound = "SMR000031368",
            concentration = "100uM",
            time = "96h",
            cellline = "HepG2",
            similarity = -0.58901,
            pValue = 0.034567890123456789012
        )
    )

    body <- list(
        status = list(
            sessionID = "test_session_12345",
            gpgene = "NA",
            gpprobe = "NA",
            Remark = "Done",
            NoOfGenes = 978L,
            NoOfProbes = "NA",
            concordanceTable = concordantsData,
            corTablePath = "/tmp/test_correlation_table.xls"
        )
    )

    httr2::response_json(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze",
        status_code = 200L,
        body = body
    )
}

# Example empty response for testing
emptyResponse <- function() {
    body <- list(
        status = list(
            sessionID = "test_session_empty",
            gpgene = "NA",
            gpprobe = "NA",
            Remark = "Done",
            NoOfGenes = 978L,
            NoOfProbes = "NA",
            concordanceTable = list("NA"),
            corTablePath = "/tmp/empty_correlation_table.xls"
        )
    )

    httr2::response_json(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze",
        status_code = 200L,
        body = body
    )
}

# Error response functions
errorResponse400 <- function() {
    httr2::response_json(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze",
        status_code = 400L,
        body = list(error = "Bad Request")
    )
}

errorResponse500 <- function() {
    httr2::response_json(
        url = "https://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze",
        status_code = 500L,
        body = list(error = "Internal Server Error")
    )
}

# Standard signature ID helpers for consistent testing
lincsKdId <- function() "LINCSKD_3935"
lincsOeId <- function() "LINCSOE_21550"
lincsCpId <- function() "LINCSCP_53465"

# VCR-compatible mock wrapper
with_mocked_responses <- function(mock_response, code) {
    # This is a placeholder for VCR integration
    # In actual usage, this would integrate with httr2's mocking system
    eval(substitute(code), envir = parent.frame())
}

# Edge case signature generators for testing
createZeroSignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = "GENE1",
        Value_LogDiffExp = 0.0
    )
}

createSmallPositiveSignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = c("GENE1", "GENE2"),
        Value_LogDiffExp = c(0.001, 0.002)
    )
}

createSmallNegativeSignature <- function() {
    tibble::tibble(
        Name_GeneSymbol = c("GENE1", "GENE2"),
        Value_LogDiffExp = c(-0.001, -0.002)
    )
}
