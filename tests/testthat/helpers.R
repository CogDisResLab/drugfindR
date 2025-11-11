################################################################################
# Simplified helper utilities for tests
# Focus: prepareSignature tests only; no file fixtures; small, deterministic
# in-memory DGE-like data.
################################################################################

## Basic seed management -------------------------------------------------------
.testSeed <- 12345L
setTestSeed <- function(seed = .testSeed) { # nocov start
    set.seed(seed)
    invisible(seed)
} # nocov end

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
        Name_GeneSymbol = genes,
        Value_LogDiffExp = stats::rnorm(nGenes, mean = 0L, sd = 1.5),
        Significance_pvalue = stats::runif(nGenes, min = 1e-4, max = 5e-2)
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

## Public test fixture accessor (signature -> DGE-like for these tests) -------
getTestFixture <- function(type = c("prepared_signature", "input_signature"), seed = NULL, ...) {
    type <- match.arg(type)
    switch(type,
        input_signature = createInputDge(seed = seed, ...),
        prepared_signature = createPreparedDge(seed = seed, ...),
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
################################################################################
