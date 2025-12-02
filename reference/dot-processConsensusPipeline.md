# Process concordants data through the complete consensus pipeline

This internal function orchestrates the complete processing pipeline for
consensus concordants analysis.

## Usage

``` r
.processConsensusPipeline(concordants, cutoff, cellLine)
```

## Arguments

- concordants:

  A combined dataframe containing all concordants data.

- cutoff:

  Numeric similarity cutoff value.

- cellLine:

  Character vector of cell lines to include, or NULL.

## Value

A processed dataframe with consensus concordants results.

## Details

This function coordinates the following processing steps:

1.  Cell line filtering (if specified)

2.  Similarity cutoff application

3.  Target grouping and maximum similarity selection

4.  Column selection and ordering

5.  Target column renaming

## Examples

``` r
if (FALSE) { # \dontrun{
testData <- data.frame(
    similarity = c(0.5, -0.8, 0.2),
    compound = c("A", "B", "C"),
    cellline = c("A375", "PC3", "A375")
)
processed <- .processConsensusPipeline(testData, 0.3, "A375")
} # }
```
