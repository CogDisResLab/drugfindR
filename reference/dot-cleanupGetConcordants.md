# Clean up temporary signature file

This internal function removes the temporary signature file created
during the getConcordants operation to prevent accumulation of temporary
files.

## Usage

``` r
.cleanupGetConcordants(signatureFile)
```

## Arguments

- signatureFile:

  Character string path to the temporary signature file to be removed.

## Value

Invisible NULL. The function is called for its side effect of removing
the temporary file.

## Details

The function checks if the specified file exists and removes it using
[`unlink()`](https://rdrr.io/r/base/unlink.html). This cleanup is
performed automatically at the end of the getConcordants operation.

## Examples

``` r
NULL
#> NULL
```
