
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drugfindR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Continuous Integration / R
Workflows](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml)
[![latest-version](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fcogdisreslab.r-universe.dev%2Fapi%2Fpackages%2FdrugfindR&query=%24.Version&style=flat&label=latest-release&color=orange)](https://github.com/CogDisResLab/drugfindR/releases/latest)
[![license](https://img.shields.io/github/license/CogDisResLab/drugfindR)](https://github.com/CogDisResLab/drugfindR/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/338354715.svg)](https://zenodo.org/badge/latestdoi/338354715)
[![Codecov test
coverage](https://codecov.io/gh/CogDisResLab/drugfindR/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/CogDisResLab/drugfindR?branch=devel)
<!-- badges: end -->

drugfindR provides R-based programmatic access to the
[iLINCS](http://www.ilincs.org) (Integrative LINCS) database for drug
repurposing and functional genomics research. The package enables
systematic analysis of gene expression signatures against comprehensive
databases of:

- **Gene Knockdown (KD)** signatures
- **Gene Overexpression (OE)** signatures
- **Chemical Perturbagen (CP)** signatures

## Overview

drugfindR streamlines transcriptomic signature analysis by:

- Querying LINCS signatures without web platform dependencies
- Processing custom transcriptomic data for concordance analysis
- Identifying candidate repurposable drugs based on signature similarity
- Discovering functional relationships between genes and compounds
- Enabling high-throughput batch processing of multiple signatures

The package integrates with standard differential expression analysis
outputs from tools like
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
and
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

## Installation

Install the stable release from
[r-universe](https://cogdisreslab.r-universe.dev/drugfindR):

``` r
install.packages("drugfindR",
    repos = c(
        "https://cogdisreslab.r-universe.dev",
        "https://cran.r-project.org"
    )
)
```

Or install the development version from
[GitHub](https://github.com/CogDisResLab/drugfindR):

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("CogDisResLab/drugfindR")
```

## Quick Start

drugfindR offers two analysis approaches:

### 1. High-Level Convenience Functions

For rapid analysis with sensible defaults, use the wrapper functions:

``` r
library(drugfindR)

# Investigate a transcriptomic signature for candidate drugs
results <- investigateSignature(
    expr = my_diffexp_data,
    outputLib = "CP", # Query chemical perturbagens
    filterThreshold = 1.5, # Filter by absolute log fold change
    geneColumn = "gene_symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)

# Investigate a specific gene target
target_results <- investigateTarget(
    target = "TP53",
    inputLib = "KD", # Use knockdown signatures
    outputLib = "CP", # Find related chemical perturbagens
    filterThreshold = 0.5
)
```

### 2. Modular Pipeline Functions

For customized workflows and fine-grained control:

``` r
# Step 1: Prepare your signature
signature <- prepareSignature(diffexp_data,
    geneColumn = "gene",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)

# Step 2: Filter by thresholds
filtered_up <- filterSignature(signature, direction = "up", threshold = 1.5)
filtered_dn <- filterSignature(signature, direction = "down", threshold = 1.5)

# Step 3: Query concordant signatures
concordants_up <- getConcordants(filtered_up, ilincsLibrary = "CP")
concordants_dn <- getConcordants(filtered_dn, ilincsLibrary = "CP")

# Step 4: Generate consensus results
consensus <- consensusConcordants(concordants_up, concordants_dn,
    paired = TRUE,
    cutoff = 0.2
)
```

## Primary Use Cases

### Drug Repurposing from Transcriptomic Data

Identify candidate drugs that reverse or mimic a disease signature by
comparing your expression data against LINCS chemical perturbagen
signatures.

### Target Validation and Functional Analysis

Discover genes and compounds with similar or opposite functional effects
to a gene of interest using knockdown and overexpression signatures.

## Documentation

- **[Full Documentation
  Website](https://cogdisreslab.github.io/drugfindR/)**
- **[Function
  Reference](https://cogdisreslab.github.io/drugfindR/reference/index.html)** -
  Complete API documentation

### Vignettes

- **[Getting
  Started](https://cogdisreslab.github.io/drugfindR/articles/getting-started.html)** -
  Quick start guide with common workflows
- **[Drug Repurposing from Transcriptomic
  Data](https://cogdisreslab.github.io/drugfindR/articles/drug-repurposing.html)** -
  Comprehensive guide for therapeutic discovery
- **[Target Investigation and Functional
  Genomics](https://cogdisreslab.github.io/drugfindR/articles/target-investigation.html)** -
  Gene and drug target analysis

### Key Functions

| Function | Purpose | Reference |
|----|----|----|
| [`investigateSignature()`](https://cogdisreslab.github.io/drugfindR/reference/investigateSignature.html) | All-in-one signature analysis | [docs](https://cogdisreslab.github.io/drugfindR/reference/investigateSignature.html) |
| [`investigateTarget()`](https://cogdisreslab.github.io/drugfindR/reference/investigateTarget.html) | Analyze a specific gene target | [docs](https://cogdisreslab.github.io/drugfindR/reference/investigateTarget.html) |
| [`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.html) | Format signatures for iLINCS | [docs](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.html) |
| [`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.html) | Apply thresholds to signatures | [docs](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.html) |
| [`getConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/getConcordants.html) | Query iLINCS for matches | [docs](https://cogdisreslab.github.io/drugfindR/reference/getConcordants.html) |
| [`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.html) | Generate consensus results | [docs](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.html) |

## Citation

If you use drugfindR in your research, please cite:

    Imami AS, Sahay S, Creeden JF, McCullumsmith RE (2024).
    drugfindR: Investigate iLINCS for candidate repurposable drugs.
    R package version 0.99.1302.
    doi:10.5281/zenodo.338354715

## Getting Help

- **Questions & Discussions**: Use [GitHub
  Discussions](https://github.com/CogDisResLab/drugfindR/discussions)
- **Bug Reports**: Open an
  [issue](https://github.com/CogDisResLab/drugfindR/issues)
- **Feature Requests**: Submit via [GitHub
  Issues](https://github.com/CogDisResLab/drugfindR/issues)

## Related Resources

- [iLINCS Portal](http://www.ilincs.org) - Web-based platform for LINCS
  data
- [LINCS Project](https://lincsproject.org/) - Library of Integrated
  Network-Based Cellular Signatures
- [O’Donovan et
  al. 2021](https://www.nature.com/articles/s41598-021-84044-9) -
  Example application for COVID-19 drug repurposing

## License

GPL-3 + file LICENSE
