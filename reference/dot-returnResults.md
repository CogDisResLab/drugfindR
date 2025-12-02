# Return results in appropriate format based on input type

This internal function formats the output results to match the input
signature type, ensuring consistent data type handling.

## Usage

``` r
.returnResults(result, inputClass)
```

## Arguments

- result:

  A tibble containing the processed results from iLINCS API.

- inputClass:

  A character vector containing the class of the original input
  signature (from `.validateGetConcordantsInput`).

## Value

The results in the appropriate format: \* S4Vectors::DataFrame if input
was a DataFrame \* tibble otherwise (for data.frame, tibble inputs)

## Details

This function ensures that the output format matches the input format
for consistency. If the original signature was provided as an
S4Vectors::DataFrame, the results are converted back to DataFrame.
Otherwise, results are returned as a tibble.

## Examples

``` r
NULL
#> NULL
```
