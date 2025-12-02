# Calculate thresholds from two threshold values

This internal function handles asymmetric filtering thresholds when two
threshold values are provided. The first value is used as the
down-regulated threshold and the second value is used as the
up-regulated threshold.

## Usage

``` r
.calculateDoubleThreshold(threshold)
```

## Arguments

- threshold:

  A numeric vector of length 2 containing the threshold values. The
  first element is the down-regulated threshold (typically negative),
  and the second element is the up-regulated threshold (typically
  positive).

## Value

A named list with two elements:

- `downThreshold`: The down-regulated threshold (`threshold[1]`)

- `upThreshold`: The up-regulated threshold (`threshold[2]`)

## Details

This function enables asymmetric filtering where different absolute
thresholds can be applied to up-regulated and down-regulated genes. This
is useful when you want to apply stricter criteria to one direction of
regulation than the other.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create asymmetric thresholds
thresholds <- .calculateDoubleThreshold(c(-2.0, 1.5))
# Returns: list(downThreshold = -2.0, upThreshold = 1.5)

# Stricter threshold for down-regulation
thresholds <- .calculateDoubleThreshold(c(-1.0, 0.5))
# Returns: list(downThreshold = -1.0, upThreshold = 0.5)

# Equal but explicit thresholds
thresholds <- .calculateDoubleThreshold(c(-1.5, 1.5))
# Returns: list(downThreshold = -1.5, upThreshold = 1.5)
} # }
```
