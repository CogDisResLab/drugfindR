# Validate getSignature input parameters

This internal function validates all input parameters for the
getSignature function to ensure they meet the required constraints.

## Usage

``` r
.validateGetSignatureInput(sigId)
```

## Arguments

- sigId:

  A character string containing the iLINCS signature ID to retrieve.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

This function performs the following validations:

- Ensures `sigId` is a character vector of length 1

- Ensures `sigId` is not empty or whitespace-only

- Validates that the signature exists in the metadata tables

## Examples

``` r
NULL
#> NULL
```
