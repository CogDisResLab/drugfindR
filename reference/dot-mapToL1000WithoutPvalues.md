# Map filtered data to L1000 format without p-values

This internal function maps the filtered expression data to the
standardized L1000 signature format, without p-value information.

## Usage

``` r
.mapToL1000WithoutPvalues(filteredData, geneColumn, logfcColumn)
```

## Arguments

- filteredData:

  A dataframe containing filtered differential expression data.

- geneColumn:

  Character string specifying the column name containing gene symbols.

- logfcColumn:

  Character string specifying the column name containing log fold-change
  values.

## Value

A tibble with standardized L1000 signature format without p-values.
