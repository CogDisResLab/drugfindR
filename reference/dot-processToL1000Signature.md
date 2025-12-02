# Process differential expression data into L1000 signature format

This internal function orchestrates the conversion of filtered
differential expression data into the standardized L1000 signature
format.

## Usage

``` r
.processToL1000Signature(
  filteredData,
  geneColumn,
  logfcColumn,
  pvalColumn = NA
)
```

## Arguments

- filteredData:

  A dataframe containing filtered differential expression data.

- geneColumn:

  Character string specifying the column name containing gene symbols.

- logfcColumn:

  Character string specifying the column name containing log fold-change
  values.

- pvalColumn:

  Character string specifying the column name containing p-values, or
  NA.

## Value

A tibble with the standardized L1000 signature format.

## Details

This function dispatches to appropriate mapping functions based on
whether p-value information is available:

1.  `.mapToL1000WithPvalues` when p-value column is specified

2.  `.mapToL1000WithoutPvalues` when p-value column is NA
