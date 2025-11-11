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

## Minimal DGE generator used by prepareSignature tests -----------------------
createDummyDge <- function(nGenes = 6L, seed = NULL) {
    if (!is.null(seed)) setTestSeed(seed)
    nGenes <- max(5L, min(10L, as.integer(nGenes)))

    # Prefer sampling from package dataset `l1000` to satisfy membership checks
    pool <- tryCatch(
        {
            if (exists("l1000", inherits = TRUE) && !is.null(l1000[["L1000"]])) unique(l1000[["L1000"]]) else NULL
        },
        error = function(e) NULL
    )
    if (is.null(pool) || length(pool) < nGenes) {
        # Fallback set of common symbols
        pool <- c("TP53", "MYC", "GAPDH", "ACTB", "EGFR", "BRCA1", "AKT1", "MAPK1", "JUN", "FOS")
    }
    genes <- sample(pool, nGenes, replace = FALSE)

    tibble::tibble(
        Name_GeneSymbol = genes,
        Value_LogDiffExp = stats::rnorm(nGenes, mean = 0, sd = 1.5),
        Significance_pvalue = stats::runif(nGenes, min = 1e-4, max = 5e-2)
    )
}

## Public test fixture accessor (signature -> DGE-like for these tests) -------
getTestFixture <- function(type = c("signature"), seed = NULL, ...) {
    type <- match.arg(type)
    switch(type,
        signature = createDummyDge(seed = seed),
        stop("Unsupported fixture type: ", type)
    )
}

################################################################################
# (Optional) Small edge-case helpers retained (lightweight) --------------------
createZeroSignature <- function() tibble::tibble(Name_GeneSymbol = "GENE_ZERO", Value_LogDiffExp = 0)
createSmallPositiveSignature <- function() tibble::tibble(Name_GeneSymbol = c("GENE_A", "GENE_B"), Value_LogDiffExp = c(0.01, 0.02))
createSmallNegativeSignature <- function() tibble::tibble(Name_GeneSymbol = c("GENE_A", "GENE_B"), Value_LogDiffExp = c(-0.01, -0.02))
################################################################################
