# Check if a single library is valid

This internal function validates whether a single library name is one of
the supported iLINCS library types.

## Usage

``` r
.validateLibrary(lib)
```

## Arguments

- lib:

  A character string containing a single library name to validate.

## Value

A logical value: TRUE if the library is valid, FALSE otherwise.

## Details

Valid library names are:

- `"CP"`: Chemical Perturbagen library

- `"KD"`: Knockdown library

- `"OE"`: Overexpression library

## Examples

``` r
NULL
#> NULL
```
