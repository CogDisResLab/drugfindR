# Apply similarity cutoff filter to concordants data

This internal function filters concordants data based on absolute
similarity values meeting or exceeding the specified cutoff threshold.

## Usage

``` r
.applySimilarityCutoff(concordants, cutoff)
```

## Arguments

- concordants:

  A dataframe containing concordants data.

- cutoff:

  Numeric similarity cutoff value.

## Value

A filtered dataframe containing only entries meeting the similarity
cutoff.

## Details

This function:

1.  Filters based on absolute similarity values

2.  Retains both positive and negative similarities above threshold

3.  Removes entries below the cutoff threshold

## Examples

``` r
if (FALSE) { # \dontrun{
testData <- data.frame(similarity = c(0.5, -0.8, 0.2, -0.1))
filtered <- .applySimilarityCutoff(testData, 0.3)
# Returns entries with |similarity| >= 0.3
} # }
```
