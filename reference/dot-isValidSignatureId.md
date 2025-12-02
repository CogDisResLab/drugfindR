# Check if a signature ID exists in the metadata tables

This internal function validates whether a signature ID exists in any of
the metadata tables (CP, KD, or OE).

## Usage

``` r
.isValidSignatureId(sigId)
```

## Arguments

- sigId:

  A character string or vector containing the signature ID(s) to
  validate.

## Value

A logical value or vector: TRUE if the signature exists, FALSE
otherwise. Returns a vector of the same length as the input when given a
vector.

## Details

This function searches all three metadata tables:

- Chemical Perturbagen (CP) metadata

- Knockdown (KD) metadata

- Overexpression (OE) metadata

The function checks the "SourceSignature" column in each metadata table
for the provided signature ID(s).

## Examples

``` r
NULL
#> NULL
```
