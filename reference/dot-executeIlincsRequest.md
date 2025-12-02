# Execute iLINCS API request with error handling

This internal function safely executes an httr2 request and captures
errors instead of raising them, allowing downstream functions to handle
them appropriately.

## Usage

``` r
.executeIlincsRequest(request, verbose = FALSE)
```

## Arguments

- request:

  An httr2_request object to be executed.

- verbose:

  Logical indicating whether to display request details. Default is
  FALSE.

## Value

An httr2_response object, including error responses that would normally
cause httr2 to raise an error.

## Details

This function configures the request to not raise errors automatically
on HTTP error status codes (4xx, 5xx) by using
[`httr2::req_error()`](https://httr2.r-lib.org/reference/req_error.html).
Instead, error responses are returned as response objects that can be
processed by
[`.processIlincsResponse()`](https://cogdisreslab.github.io/drugfindR/reference/dot-processIlincsResponse.md)
to generate appropriate error messages with context.

The function handles:

- Network connection errors

- HTTP error status codes (400, 401, 403, 404, 500, etc.)

- Timeout errors

- Other httr2 request failures

## Examples

``` r
NULL
#> NULL
```
