# Validate getConcordants input parameters

This internal function validates the input parameters for the
getConcordants function to ensure they meet the required constraints.

## Usage

``` r
.validateGetConcordantsInput(signature, ilincsLibrary)
```

## Arguments

- signature:

  A data.frame-like object containing the signature data. Must have the
  required iLINCS signature structure with columns: signatureID,
  ID_geneid, Name_GeneSymbol, Value_LogDiffExp, Significance_pvalue.

- ilincsLibrary:

  Character string specifying the iLINCS library to search. Must be one
  of "OE", "KD", or "CP".

## Value

A character vector containing the class of the input signature. This is
used internally to determine the return type format.

## Details

This function performs comprehensive validation:

1.  Ensures `signature` is a data.frame-like object (data.frame, tibble,
    or S4Vectors::DataFrame)

2.  Validates complete signature structure via
    [`stopIfInvalidSignature()`](https://cogdisreslab.github.io/drugfindR/reference/stopIfInvalidSignature.md)

3.  Validates `ilincsLibrary` is one of the supported libraries

The signature must conform to the iLINCS expected structure. Use
[`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.md)
to ensure proper formatting.

## See also

`[ stopIfInvalidSignature() ]` for signature structure validation,
`[ prepareSignature() ]` for signature preparation

## Examples

``` r
NULL
#> NULL
```
