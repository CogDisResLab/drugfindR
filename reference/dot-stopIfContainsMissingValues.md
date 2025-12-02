# Validate signature for missing values

This internal function checks if the signature data frame contains any
missing (NA) values, which are not allowed in iLINCS signature data.

## Usage

``` r
.stopIfContainsMissingValues(signature)
```

## Arguments

- signature:

  A data.frame-like object containing signature data.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

The function scans the entire signature data frame for missing values.
iLINCS requires complete data for all signature analysis, so any NA
values will cause the function to stop with an informative error message
indicating which columns contain missing values.

## Examples

``` r
NULL
#> NULL
```
