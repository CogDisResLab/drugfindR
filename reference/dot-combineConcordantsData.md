# Combine concordants dataframes for consensus analysis

This internal function combines one or more concordants dataframes into
a single dataframe for further processing.

## Usage

``` r
.combineConcordantsData(dots)
```

## Arguments

- dots:

  A list of dataframes to combine.

## Value

A combined dataframe with all input data.

## Details

This function:

1.  Combines multiple dataframes using row binding

2.  Preserves all columns from input dataframes

3.  Handles cases where dataframes have different column sets

## Examples

``` r
if (FALSE) { # \dontrun{
df1 <- data.frame(similarity = 0.5, compound = "A")
df2 <- data.frame(similarity = -0.3, compound = "B")
combined <- .combineConcordantsData(list(df1, df2))
} # }
```
