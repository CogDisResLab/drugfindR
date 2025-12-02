# Return the internal iLINCS Library ID for a given library

This internal function maps user-friendly library names to the internal
library identifiers used by the iLINCS API.

## Usage

``` r
.returnLibrary(lib)
```

## Arguments

- lib:

  A character string specifying the library name. Must be one of "OE",
  "KD", or "CP".

## Value

A character string containing the corresponding iLINCS library ID.

## Details

The mapping between user library names and iLINCS internal IDs is:

- `"OE"` (Overexpression) -\> `"LIB_11"`

- `"KD"` (Knockdown) -\> `"LIB_6"`

- `"CP"` (Chemical Perturbagen) -\> `"LIB_5"`

The function validates the input library name before mapping and will
stop execution if an invalid library is provided.

## See also

`[ stopIfInvalidLibraries() ]` for library validation details

## Examples

``` r
NULL
#> NULL
```
