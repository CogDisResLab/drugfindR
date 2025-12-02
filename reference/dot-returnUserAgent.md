# Return a string suitable as a User-Agent for the iLINCS API

This internal function constructs a standardized User-Agent string for
HTTP requests to the iLINCS API, including package name, version, and
repository URL for identification and debugging purposes.

## Usage

``` r
.returnUserAgent()
```

## Value

A character string formatted as a User-Agent header value.

## Details

The User-Agent string follows the format:
`"drugfindR/<current version>; https://github.com/CogDisResLab/drugfindR"`

This helps iLINCS administrators identify requests from this package and
assists with debugging if issues arise. The version is automatically
retrieved from the package metadata.

## Examples

``` r
NULL
#> NULL
```
