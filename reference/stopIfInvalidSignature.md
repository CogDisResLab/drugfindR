# Validate signature data structure and content

This function performs comprehensive validation of signature data to
ensure it meets the requirements for iLINCS analysis.

## Usage

``` r
stopIfInvalidSignature(signature)
```

## Arguments

- signature:

  A data.frame-like object containing signature data that should be
  validated for iLINCS compatibility.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

This function performs two main validation checks:

1.  Column structure validation via
    [`.stopIfInvalidColNames()`](https://cogdisreslab.github.io/drugfindR/reference/dot-stopIfInvalidColNames.md)

2.  Missing value validation via `[.stopIfContainsMissingValues()]`

The signature must have exactly the required columns in the correct
order and cannot contain any missing (NA) values.

## See also

`[ prepareSignature() ]` for preparing signatures that meet these
requirements, `[ .stopIfInvalidColNames() ]` for column validation
details, `[ .stopIfContainsMissingValues() ]` for missing value
validation details

## Examples

``` r
NULL
#> NULL
```
