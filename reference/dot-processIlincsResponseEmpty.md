# Handle empty concordance table responses

This internal function creates an empty tibble with the correct
structure when the iLINCS API returns no concordant signatures.

## Usage

``` r
.processIlincsResponseEmpty(sigDirection, ilincsLibrary)
```

## Arguments

- sigDirection:

  Character string indicating the signature direction.

- ilincsLibrary:

  Character string specifying the iLINCS library used.

## Value

A tibble with zero rows and the correct column structure.
