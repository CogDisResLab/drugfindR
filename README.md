<!-- README.md is generated from README.Rmd. Please edit that file -->

# drugfindR: An R package to search iLINCS databases for small molecules [![Lifecycle: stable](man/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: start -->

[![drugfindR status badge](https://cogdisreslab.r-universe.dev/badges/drugfindR)](https://cogdisreslab.r-universe.dev/drugfindR)
![GitHub R package version](https://img.shields.io/github/r-package/v/CogDisResLab/drugfindR)
[![latest-version](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fcogdisreslab.r-universe.dev%2Fapi%2Fpackages%2FdrugfindR&query=%24.Version&style=flat&label=latest-release&color=orange)](https://github.com/CogDisResLab/drugfindR/releases/latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8346520.svg)](https://doi.org/10.5281/zenodo.8346520)
[![license](https://img.shields.io/github/license/CogDisResLab/drugfindR)](https://github.com/CogDisResLab/drugfindR/blob/main/LICENSE)

[![Continuous Integration / R Workflows](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/CogDisResLab/drugfindR/actions/workflows/rworkflows.yml)
[![Codecov test coverage](https://codecov.io/gh/CogDisResLab/drugfindR/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/CogDisResLab/drugfindR?branch=devel)

![Libraries.io dependency status for GitHub repo](https://img.shields.io/librariesio/github/CogDisResLab/drugfindR)
[![CodeFactor](https://www.codefactor.io/repository/github/cogdisreslab/drugfindr/badge)](https://www.codefactor.io/repository/github/cogdisreslab/drugfindr)
![REUSE Compliance](https://img.shields.io/reuse/compliance/github.com%2FCogDisResLab%2FdrugfindR.git)

![GitHub repo size](https://img.shields.io/github/repo-size/CogDisResLab/drugfindR)
![GitHub Issues](https://img.shields.io/github/issues/CogDisResLab/drugfindR)
![GitHub Pull Requests](https://img.shields.io/github/issues-pr/CogDisResLab/drugfindR)

<!-- badges: end -->

drugfindR allows convenient access to the iLINCS Gene Knockdown, Gene
Overexpression and Chemical Perturbagen databases and allows you to
generate and investigate signatures to identify relevant genes and
drugs.

## Installation

You can install the released version of drugfindR from
[r-universe](https://cogdisreslab.r-universe.dev/drugfindR) with:

```r
install.packages("drugfindR",
    repos = c(
        "https://cogdisreslab.r-universe.dev",
        "https://cran.r-project.org"
    )
)
```

And the development version from [GitHub](https://github.com/) with:

```r
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
