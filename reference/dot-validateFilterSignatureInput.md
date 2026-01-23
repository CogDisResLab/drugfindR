# Validate filterSignature input parameters

This internal function validates all input parameters for the
filterSignature function to ensure they meet the required constraints
and are mutually compatible.

## Usage

``` r
.validateFilterSignatureInput(signature, direction, threshold, prop)
```

## Arguments

- signature:

  A data.frame-like object (data.frame, tibble, or DataFrame) containing
  the L1000 signature data.

- direction:

  Character string specifying the filtering direction. Must be one of
  "up", "down", or "any".

- threshold:

  Numeric value or vector specifying absolute threshold(s). Can be NULL,
  a single value, or a vector of two values. Cannot be specified
  together with `prop`.

- prop:

  Numeric value specifying the proportion for quantile-based filtering.
  Must be between 0 and 1. Cannot be specified together with
  `threshold`.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

This function performs the following validations in order:

1.  Ensures `signature` is a data.frame-like object

2.  Validates `direction` is one of the allowed values

3.  Verifies that only one of `threshold` or `prop` is specified

4.  For `threshold`: checks length (1-2 values) and order (lower,
    higher)

5.  For `prop`: checks it's a single value between 0 and 1

## Examples

``` r
if (FALSE) { # \dontrun{
# Valid calls (no errors)
sig <- data.frame(Value_LogDiffExp = c(-2, -1, 0, 1, 2))
.validateFilterSignatureInput(sig, "any", 1.0, NULL)
.validateFilterSignatureInput(sig, "up", NULL, 0.1)
.validateFilterSignatureInput(sig, "down", c(-1.5, 1.0), NULL)

# Invalid calls (will throw errors)
.validateFilterSignatureInput(sig, "invalid", 1.0, NULL) # Invalid direction
.validateFilterSignatureInput(sig, "any", 1.0, 0.1) # Both threshold and prop
.validateFilterSignatureInput(sig, "any", NULL, NULL) # Neither threshold nor prop
.validateFilterSignatureInput(sig, "any", c(1, 2, 3), NULL) # Too many thresholds
.validateFilterSignatureInput(sig, "any", c(2, 1), NULL) # Wrong threshold order
.validateFilterSignatureInput(sig, "any", NULL, 1.5) # Proportion > 1
.validateFilterSignatureInput(sig, "any", NULL, -0.1) # Proportion < 0
} # }
```
