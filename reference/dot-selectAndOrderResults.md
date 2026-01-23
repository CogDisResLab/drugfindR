# Select and order consensus results columns

This internal function selects the relevant columns for consensus
results and orders them appropriately for output.

## Usage

``` r
.selectAndOrderResults(concordants)
```

## Arguments

- concordants:

  A dataframe containing processed concordants data.

## Value

A dataframe with selected and ordered columns for consensus output.

## Details

This function:

1.  Selects standard consensus output columns

2.  Orders results by descending absolute similarity

3.  Handles both CP/KD libraries (with concentration) and OE libraries
    (without)

## Examples

``` r
if (FALSE) { # \dontrun{
testData <- data.frame(
    signatureid = "SIG1",
    compound = "A",
    cellline = "A375",
    similarity = 0.8,
    sig_direction = "Up",
    pValue = 0.01
)
selected <- .selectAndOrderResults(testData)
} # }
```
