# Handle API error responses for signature retrieval

This internal function processes error responses from the iLINCS API and
generates appropriate error messages for signature retrieval failures.

## Usage

``` r
.processSignatureResponseError(response)
```

## Arguments

- response:

  An httr2 response object from the iLINCS API.

## Value

This function always stops execution with an error message.
