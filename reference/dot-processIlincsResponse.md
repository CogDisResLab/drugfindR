# Process iLINCS API response into concordant signatures

This internal function dispatches to appropriate handlers based on the
response status and content from the iLINCS API.

## Usage

``` r
.processIlincsResponse(response, sigDirection, ilincsLibrary)
```

## Arguments

- response:

  An httr2 response object from the iLINCS API.

- sigDirection:

  Character string indicating the signature direction ("Up", "Down", or
  "Any").

- ilincsLibrary:

  Character string specifying the iLINCS library used ("CP", "KD", or
  "OE").

## Value

A tibble containing concordant signature data with standardized column
names and rounded numerical values.

## Details

The function dispatches to specialized handlers:

1.  `.processIlincsResponseError` for HTTP error responses

2.  `.processIlincsResponseEmpty` for empty concordance tables

3.  `.processIlincsResponseSuccess` for successful responses with data

The resulting tibble always contains these columns in order:

- `signatureid`: Unique signature identifier

- `treatment`: Drug/treatment name (compound renamed for CP library)

- `concentration`: Drug concentration (NA for KD/OE libraries)

- `time`: Treatment duration

- `cellline`: Cell line used

- `similarity`: Similarity score (rounded to 8 decimal places)

- `pValue`: Statistical significance (rounded to 20 decimal places)

- `sig_direction`: Signature direction ("Up", "Down", or "Any")

- `sig_type`: Library type description
