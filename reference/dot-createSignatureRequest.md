# Create HTTP request for iLINCS signature retrieval

This internal function constructs and configures the HTTP request object
for retrieving signature data from the iLINCS API.

## Usage

``` r
.createSignatureRequest(sigId)
```

## Arguments

- sigId:

  A character string containing the iLINCS signature ID to retrieve.

## Value

An httr2 request object configured for the iLINCS downloadSignature
endpoint.

## Details

This function builds a complete HTTP request by:

1.  Setting the base URL using
    [`.ilincsBaseUrl()`](https://cogdisreslab.github.io/drugfindR/reference/dot-ilincsBaseUrl.md)

2.  Appending the API path: "ilincsR/downloadSignature"

3.  Adding query parameters: sigID and noOfTopGenes (set to Inf for all
    genes)

4.  Setting the HTTP method to POST

5.  Adding a user agent string using
    [`.returnUserAgent()`](https://cogdisreslab.github.io/drugfindR/reference/dot-returnUserAgent.md)

The request is configured but not executed - it must be performed using
the request execution function.

## Examples

``` r
NULL
#> NULL
```
