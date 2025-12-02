# Get the L1000 Signature from iLINCS **\[stable\]**

This function acts as the entrypoint to the iLINCS database. This takes
in an ID and returns the signature after making a call to the iLINCS
database. The function automatically detects whether the signature is an
L1000 signature based on the signature ID and metadata tables, and
retrieves all available genes for comprehensive signature analysis.

## Usage

``` r
getSignature(sigId)
```

## Arguments

- sigId:

  character. The ilincs signature_id

## Value

a tibble with the signature data containing the following columns: \*
`signatureID`: The signature identifier \* `ID_geneid`: Gene IDs
(Entrez) \* `Name_GeneSymbol`: Gene symbols \* `Value_LogDiffExp`: Log
fold-change values \* `Significance_pvalue`: Statistical significance
p-values

## Examples

``` r
# Input validation example (no API call)
# Demonstrates proper signature ID format validation
tryCatch(
    getSignature(""), # Empty string should error
    error = function(e) message("Expected error: empty signature ID")
)
#> Expected error: empty signature ID

# \donttest{
# These examples require network access to the iLINCS API

# Get the L1000 signature for LINCSKD_28
kdSignature <- getSignature("LINCSKD_28")
head(kdSignature)
#> # A tibble: 6 × 5
#>   signatureID ID_geneid Name_GeneSymbol Value_LogDiffExp Significance_pvalue
#>   <chr>       <chr>     <chr>                      <dbl>               <dbl>
#> 1 LINCSKD_28  6714      SRC                         5.06            0.000126
#> 2 LINCSKD_28  23039     XPO7                       -4.84            0.000857
#> 3 LINCSKD_28  5300      PIN1                       -4.47            0.00131 
#> 4 LINCSKD_28  51335     NGRN                       -3.97            0.00169 
#> 5 LINCSKD_28  1465      CSRP1                      -4.97            0.00300 
#> 6 LINCSKD_28  23149     FCHO1                       3.97            0.00378 

# Get an overexpression signature (L1000 status is automatically detected)
oeSignature <- getSignature("LINCSOE_1000")
head(oeSignature)
#> # A tibble: 6 × 5
#>   signatureID  ID_geneid Name_GeneSymbol Value_LogDiffExp Significance_pvalue
#>   <chr>        <chr>     <chr>                      <dbl>               <dbl>
#> 1 LINCSOE_1000 54733     SLC35F2                    1.67              0.00110
#> 2 LINCSOE_1000 9961      MVP                       -2.33              0.00126
#> 3 LINCSOE_1000 1647      GADD45A                    1.14              0.00167
#> 4 LINCSOE_1000 5985      RFC5                       2.71              0.00273
#> 5 LINCSOE_1000 8508      NIPSNAP1                  -1.11              0.00308
#> 6 LINCSOE_1000 54541     DDIT4                     -0.985             0.00341

# Check the structure of retrieved signature
str(kdSignature)
#> tibble [978 × 5] (S3: tbl_df/tbl/data.frame)
#>  $ signatureID        : chr [1:978] "LINCSKD_28" "LINCSKD_28" "LINCSKD_28" "LINCSKD_28" ...
#>  $ ID_geneid          : chr [1:978] "6714" "23039" "5300" "51335" ...
#>  $ Name_GeneSymbol    : chr [1:978] "SRC" "XPO7" "PIN1" "NGRN" ...
#>  $ Value_LogDiffExp   : num [1:978] 5.06 -4.84 -4.47 -3.97 -4.97 ...
#>  $ Significance_pvalue: num [1:978] 0.000126 0.000857 0.001314 0.001691 0.003001 ...
# }
```
