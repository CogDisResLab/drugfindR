# Filter concordants data by cell line

This internal function filters the concordants data to include only the
specified cell lines.

## Usage

``` r
.filterByCellLine(concordants, cellLine)
```

## Arguments

- concordants:

  A dataframe containing concordants data.

- cellLine:

  A character vector of cell lines to include, or NULL for no filtering.

## Value

A filtered dataframe containing only the specified cell lines.

## Details

This function:

1.  Filters data based on the cellline column

2.  Returns original data if cellLine is NULL

3.  Handles cases where no data matches the specified cell lines

## Examples

``` r
if (FALSE) { # \dontrun{
testData <- data.frame(
    similarity = c(0.5, -0.3, 0.7),
    cellline = c("A375", "PC3", "MCF7")
)
filtered <- .filterByCellLine(testData, c("A375", "PC3"))
} # }
```
