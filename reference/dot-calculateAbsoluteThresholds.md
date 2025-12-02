# Calculate thresholds using absolute threshold values

This internal function coordinates the calculation of filtering
thresholds when absolute threshold values are provided. It dispatches to
the appropriate calculation function based on the number of threshold
values provided.

## Usage

``` r
.calculateAbsoluteThresholds(threshold)
```

## Arguments

- threshold:

  A numeric value or vector specifying the absolute threshold(s).

  Can be:

      - A single value: Dispatched to [.calculateSingleThreshold()]
      - A vector of two values: Dispatched to [.calculateDoubleThreshold()]

## Value

A named list with two elements:

- `downThreshold`: The threshold for down-regulated genes

- `upThreshold`: The threshold for up-regulated genes

## Details

This function serves as a dispatcher that:

- Checks the length of the threshold parameter

- Calls the appropriate threshold calculation function

- Throws an error if an invalid number of thresholds is provided

The function ensures that only single values or pairs of values are
accepted, maintaining the integrity of the filtering logic.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single threshold - creates symmetric thresholds
thresholds <- .calculateAbsoluteThresholds(1.0)
# Returns: list(downThreshold = -1.0, upThreshold = 1.0)

# Double threshold - uses provided values
thresholds <- .calculateAbsoluteThresholds(c(-1.5, 2.0))
# Returns: list(downThreshold = -1.5, upThreshold = 2.0)

# Invalid - too many values (will throw error)
# thresholds <- .calculateAbsoluteThresholds(c(1.0, 2.0, 3.0))
} # }
```
