# Apply target column renaming to consensus results

This internal function applies the standard target column renaming to
produce the final consensus concordants output format.

## Usage

``` r
.applyTargetRenaming(concordants)
```

## Arguments

- concordants:

  A dataframe containing selected consensus results.

## Value

A dataframe with renamed columns following consensus output standards.

## Details

This function:

1.  Applies targetRename function to standardize column names

2.  Converts internal column names to user-facing consensus format

3.  Handles different library types appropriately

## Examples

``` r
testData <- data.frame(
    signatureid = "SIG1",
    compound = "A",
    cellline = "A375",
    similarity = 0.8
)
renamed <- .applyTargetRenaming(testData)
#> Error in .applyTargetRenaming(testData): could not find function ".applyTargetRenaming"
```
