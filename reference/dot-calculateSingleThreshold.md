# Calculate thresholds from single threshold value

This internal function creates symmetric filtering thresholds from a
single threshold value. The input value is used as the positive
threshold, and its negative is used as the negative threshold.

## Usage

``` r
.calculateSingleThreshold(threshold)
```

## Arguments

- threshold:

  A single positive numeric value representing the absolute threshold
  for filtering.

## Value

A named list with two elements:

- `downThreshold`: The negative threshold (-threshold)

- `upThreshold`: The positive threshold (threshold)

## Details

This function is used when a single threshold value is provided to
`filterSignature`. It creates symmetric thresholds where genes with log
fold-change values greater than or equal to the positive threshold
(up-regulated) or less than or equal to the negative threshold
(down-regulated) are retained.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create symmetric thresholds from threshold = 1.5
thresholds <- .calculateSingleThreshold(1.5)
# Returns: list(downThreshold = -1.5, upThreshold = 1.5)

# Create symmetric thresholds from threshold = 0.8
thresholds <- .calculateSingleThreshold(0.8)
# Returns: list(downThreshold = -0.8, upThreshold = 0.8)
} # }
```
