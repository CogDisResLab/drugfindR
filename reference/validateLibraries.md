# Check if multiple libraries are valid

This function validates whether all provided library names are supported
iLINCS library types.

## Usage

``` r
validateLibraries(libs)
```

## Arguments

- libs:

  A character vector of library names to validate.

## Value

A logical value: TRUE if all libraries are valid, FALSE if any are
invalid.

## Details

This function uses
[`.validateLibrary()`](https://cogdisreslab.github.io/drugfindR/reference/dot-validateLibrary.md)
to check each library individually and returns TRUE only if all
libraries are valid. It's used internally to validate library parameters
before API calls.

## See also

`[ .validateLibrary() ]` for single library validation,
`[ stopIfInvalidLibraries() ]` for validation with error handling

## Examples

``` r
NULL
#> NULL
```
