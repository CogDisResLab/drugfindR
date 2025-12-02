# Validate signature column names

This internal function checks if the signature data frame has the
expected column names in the correct order for iLINCS compatibility.

## Usage

``` r
.stopIfInvalidColNames(signature)
```

## Arguments

- signature:

  A data.frame-like object containing signature data.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

The function validates that the signature has exactly the following
columns in the specified order:

1.  `signatureID`: Signature identifier

2.  `ID_geneid`: Gene ID

3.  `Name_GeneSymbol`: Gene symbol

4.  `Value_LogDiffExp`: Log fold-change expression value

5.  `Significance_pvalue`: Statistical significance p-value

## Examples

``` r
NULL
#> NULL
```
