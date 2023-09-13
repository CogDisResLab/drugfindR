
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drugfindR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![Continuous Integration / R
Workflows](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml)
[![codecov](https://codecov.io/gh/CogDisResLab/drugfindR/branch/main/graph/badge.svg?token=FeAvIeTAiz)](https://codecov.io/gh/CogDisResLab/drugfindR)
![GitHub release (latest SemVer including
pre-releases)](https://img.shields.io/github/v/release/CogDisResLab/drugfindR?include_prereleases&label=latest-release)
[![GitHub
license](https://img.shields.io/github/license/CogDisResLab/drugfindR)](https://github.com/CogDisResLab/drugfindR/blob/main/LICENSE)
<!-- badges: end -->

drugfindR allows convenient access to the iLINCS Gene Knockdown, Gene
Overexpression and Chemical Perturbagen databases and allows you to
generate and investigate signatures to identify relevant genes and
drugs.

## Installation

You can install the released version of drugfindR from
[bioconductor](https://bioconductor.org/) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("drugfindR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("CogDisResLab/drugfindR")
```
