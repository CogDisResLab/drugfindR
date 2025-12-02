# Stop if the libraries are invalid

This internal function validates library specifications and stops
execution with an informative error message if any invalid libraries are
found.

## Usage

``` r
stopIfInvalidLibraries(libs)
```

## Arguments

- libs:

  A character vector of library names to validate. Each library must be
  one of "OE", "KD", or "CP".

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

This function validates that all provided library names are supported:

- `"OE"`: Overexpression library (LIB_11)

- `"KD"`: Knockdown library (LIB_6)

- `"CP"`: Chemical Perturbagen library (LIB_5)

If any invalid libraries are found, the function provides a detailed
error message listing the invalid libraries and the expected options.

## See also

`[ validateLibraries() ]` for the underlying validation logic,
`[ .validateLibrary() ]` for single library validation

## Examples

``` r
NULL
#> NULL
```
