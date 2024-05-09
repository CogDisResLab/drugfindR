
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drugfindR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![Continuous Integration / R
Workflows](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml)
[![codecov](https://codecov.io/gh/CogDisResLab/drugfindR/branch/main/graph/badge.svg?token=FeAvIeTAiz)](https://codecov.io/gh/CogDisResLab/drugfindR)
\[![latest-version](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fcogdisreslab.r-universe.dev%2Fapi%2Fpackages%2FdrugfindR&query=%24.Version&style=flat&label=latest-release&color=orange)\]
[![license](https://img.shields.io/github/license/CogDisResLab/drugfindR)](https://github.com/CogDisResLab/drugfindR/blob/main/LICENSE)

[![DOI](https://zenodo.org/badge/338354715.svg)](https://zenodo.org/badge/latestdoi/338354715)
<!-- badges: end -->

drugfindR allows convenient access to the iLINCS Gene Knockdown, Gene
Overexpression and Chemical Perturbagen databases and allows you to
generate and investigate signatures to identify relevant genes and
drugs.

## Installation

You can install the released version of drugfindR from
[r-universe](https://cogdisreslab.r-universe.dev/drugfindR) with:

``` r
install.packages("drugfindR",
    repos = c(
        "https://cogdisreslab.r-universe.dev",
        "https://cran.r-project.org"
    )
)
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
    (`investigateSignatures` or `investigateTarget`) to generate
    results. This approach uses the building block functions under the
    hood with sensible defaults and returns a final result that can be
    used for further analysis.
