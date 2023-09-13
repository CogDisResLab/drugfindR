
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

## Usage

This package has two different interfaces that can be used to generate
the results:

1.  The first way is to use the individual building block functions to
    generate results. This is useful if you want to use the results in
    your own analysis pipeline or want more control over the results.

2.  The second way is to use one of the convenience functions
    (`investigate_signatures` or `investigate_target`) to generate
    results. This approach uses the building block functions under the
    hood with sensible defaults and returns a final result that can be
    used for further analysis.
