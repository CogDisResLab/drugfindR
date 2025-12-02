# Load the correct metadata table for a given library

This internal function retrieves the appropriate metadata table based on
the specified iLINCS library type.

## Usage

``` r
.loadMetadata(lib)
```

## Arguments

- lib:

  A character string specifying the library type. Must be one of "OE",
  "KD", or "CP".

## Value

A tibble containing the metadata for the specified library. The
structure varies by library type but typically includes columns for
signature identifiers, treatments, cell lines, and other metadata.

## Details

The function loads pre-compiled metadata tables for each library:

- `"OE"`: Overexpression metadata (`oeMetadata`)

- `"KD"`: Knockdown metadata (`kdMetadata`)

- `"CP"`: Chemical Perturbagen metadata (`cpMetadata`)

These metadata tables are included with the package and contain
information about available signatures in each iLINCS library.

## Examples

``` r
NULL
#> NULL
```
