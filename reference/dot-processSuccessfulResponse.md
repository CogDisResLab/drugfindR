# Process successful API response into signature data frame

This internal function processes a successful HTTP response from the
iLINCS API and converts it into a standardized signature data frame.

## Usage

``` r
.processSuccessfulResponse(response)
```

## Arguments

- response:

  An httr2 response object from a successful iLINCS API call.

## Value

A tibble containing the signature data with standardized columns: \*
signatureID: The signature identifier \* ID_geneid: Character gene
identifiers \* Name_GeneSymbol: Gene symbols \* Value_LogDiffExp: Log
fold-change values (rounded to 12 decimal places) \*
Significance_pvalue: P-values (rounded to 12 decimal places)

## Details

This function:

1.  Extracts JSON data from the response body

2.  Maps the "signature" elements from the response

3.  Flattens the nested structure into a data frame

4.  Removes the "PROBE" column (not needed for analysis)

5.  Converts gene IDs to character format

6.  Rounds numeric values to 12 decimal places for consistency

7.  Adds signature metadata including L1000 status

The rounding ensures consistent precision across different platforms and
prevents floating-point precision issues in downstream analyses.

## Examples

``` r
NULL
#> NULL
```
