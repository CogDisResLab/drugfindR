# Handle iLINCS API response errors

This internal function processes error responses from the iLINCS API and
generates appropriate error messages.

## Usage

``` r
.processIlincsResponseError(response)
```

## Arguments

- response:

  An httr2 response object from the iLINCS API.

## Value

This function always stops execution with an error message.
