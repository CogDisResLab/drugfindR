# Detect signature direction from expression values

This internal function analyzes the log fold-change values in a
signature to determine the overall direction of regulation.

## Usage

``` r
.detectSignatureDirection(signature)
```

## Arguments

- signature:

  A data.frame-like object containing the signature data. Must have a
  column named "Value_LogDiffExp" with log fold-change values.

## Value

Character string indicating signature direction: "Up", "Down", or "Any".

## Details

The function examines the "Value_LogDiffExp" column to determine
direction:

- "Up": All expression values are greater than or equal to zero

- "Down": All expression values are less than or equal to zero

- "Any": Mixed positive and negative values

Note that zero values are considered "Up" direction. This direction
information is used by iLINCS for signature analysis and is included in
the output results.

## Examples

``` r
NULL
#> NULL
```
