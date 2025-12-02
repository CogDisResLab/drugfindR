# Prepare signature file for iLINCS upload

This internal function creates a temporary file containing the signature
data formatted for upload to the iLINCS API.

## Usage

``` r
.prepareSignatureFile(signature)
```

## Arguments

- signature:

  A data.frame-like object containing the signature data.

## Value

Character string path to the temporary signature file.

## Details

The function creates a temporary file with a ".xls" extension and writes
the signature data as a tab-separated file. The file is automatically
cleaned up by the system when the R session ends.

## Examples

``` r
NULL
#> NULL
```
