# Handle successful concordance table responses

This internal function processes successful responses from the iLINCS
API and formats the concordant signature data with standardized columns.

## Usage

``` r
.processIlincsResponseSuccess(concordanceTables, sigDirection, ilincsLibrary)
```

## Arguments

- concordanceTables:

  List containing concordance table data from API response.

- sigDirection:

  Character string indicating the signature direction.

- ilincsLibrary:

  Character string specifying the iLINCS library used.

## Value

A tibble containing processed concordant signature data.
