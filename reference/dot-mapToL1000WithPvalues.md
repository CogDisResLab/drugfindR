# Map filtered data to L1000 format with p-values

This internal function maps the filtered differential expression data to
the standardized L1000 signature format, including p-value information.

## Usage

``` r
.mapToL1000WithPvalues(filteredData, geneColumn, logfcColumn, pvalColumn)
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

  Character string specifying the column name containing p-values.

## Value

A tibble with standardized L1000 signature format including p-values.
