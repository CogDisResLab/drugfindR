# Process iLINCS API response for signature retrieval

This internal function dispatches to appropriate handlers based on the
response status from the iLINCS API.

## Usage

``` r
.processSignatureResponse(response)
```

## Arguments

- response:

  An httr2 response object from the iLINCS API.

## Value

A tibble containing signature data with standardized columns.

## Details

The function dispatches to specialized handlers:

1.  [`.processSignatureResponseError()`](https://cogdisreslab.github.io/drugfindR/reference/dot-processSignatureResponseError.md)
    for HTTP error responses

2.  [`.processSuccessfulResponse()`](https://cogdisreslab.github.io/drugfindR/reference/dot-processSuccessfulResponse.md)
    for successful responses with data

The resulting tibble contains these columns:

- `signatureID`: The signature identifier

- `ID_geneid`: Character gene identifiers

- `Name_GeneSymbol`: Gene symbols

- `Value_LogDiffExp`: Log fold-change values (rounded to 12 decimal
  places)

- `Significance_pvalue`: P-values (rounded to 12 decimal places)
