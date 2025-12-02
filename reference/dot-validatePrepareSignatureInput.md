# Validate prepareSignature input parameters

This internal function validates all input parameters for the
prepareSignature function to ensure they meet the required constraints.

## Usage

``` r
.validatePrepareSignatureInput(dge, geneColumn, logfcColumn, pvalColumn)
```

## Arguments

- dge:

  A dataframe-like object containing differential gene expression data.

- geneColumn:

  Character string specifying the column name containing gene symbols.

- logfcColumn:

  Character string specifying the column name containing log fold-change
  values.

- pvalColumn:

  Character string specifying the column name containing p-values, or
  NA.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

This function performs the following validations:

1.  Ensures all column names are character strings

2.  Validates that specified columns exist in the input dataframe

3.  Checks that the dataframe is not empty

## Examples

``` r
NULL
#> NULL
```
