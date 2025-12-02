# Create iLINCS API request

This internal function constructs and executes the HTTP request to the
iLINCS API for concordant signature analysis.

## Usage

``` r
.generateIlincsRequest(signatureFile, ilincsLibrary)
```

## Arguments

- signatureFile:

  Character string path to the signature file to upload.

- ilincsLibrary:

  Character string specifying the iLINCS library to search. Must be one
  of "OE", "KD", or "CP".

## Value

An httr2 response object from the iLINCS API.

## Details

The function:

1.  Maps the library name to the internal iLINCS library ID

2.  Constructs a multipart POST request with the signature file

3.  Includes appropriate user agent and API endpoint

4.  Executes the request and returns the response

The library mapping is:

- CP (Chemical Perturbagen): LIB_5

- KD (Knockdown): LIB_6

- OE (Overexpression): LIB_11

## Examples

``` r
NULL
#> NULL
```
