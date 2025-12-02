# Compute consensus concordant signatures from a single input signature

This internal helper wraps the common paired / unpaired workflow used by
[`investigateSignature()`](https://cogdisreslab.github.io/drugfindR/reference/investigateSignature.md)
and
[`investigateTarget()`](https://cogdisreslab.github.io/drugfindR/reference/investigateTarget.md)
for a single already prepared or retrieved signature. It applies
directional filtering, queries iLINCS for concordant signatures, and
collapses results via
[`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md).

## Usage

``` r
.computeConsensusFromSignature(
  signature,
  outputLib,
  filterThreshold = NULL,
  filterProp = NULL,
  similarityThreshold = 0.321,
  paired = TRUE,
  outputCellLines = NULL
)
```

## Arguments

- signature:

  A data.frame / tibble / DataFrame produced by
  [`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.md)
  or
  [`getSignature()`](https://cogdisreslab.github.io/drugfindR/reference/getSignature.md)
  with standard signature columns.

- outputLib:

  Character. One of "OE", "KD", or "CP" indicating the iLINCS library to
  search for concordant signatures.

- filterThreshold:

  Numeric (optional). Absolute threshold(s) passed to
  [`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md).
  Use either `filterThreshold` or `filterProp`.

- filterProp:

  Numeric in (0, 0.5\] (optional). Proportion for quantile based
  filtering in
  [`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md).
  Ignored if `filterThreshold` is supplied.

- similarityThreshold:

  Numeric in 0..1. Minimum absolute similarity retained by
  [`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md).

- paired:

  Logical. If TRUE perform separate up / down filtering and concordance
  queries; otherwise aggregate direction = "any".

- outputCellLines:

  Optional character vector restricting target cell lines during
  consensus filtering. Passed to
  [`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md).

## Value

A tibble of consensus concordant signatures with standardized target
columns (already renamed via internal consensus pipeline). Columns
include `TargetSignature`, `Target`, `TargetCellLine`, `Similarity`,
`pValue`, `InputSigDirection`, `SignatureType`, and optional time /
concentration.

## Details

Error handling is delegated to component functions:

- Library validation via
  [`stopIfInvalidLibraries()`](https://cogdisreslab.github.io/drugfindR/reference/stopIfInvalidLibraries.md)

- Signature structure via
  [`stopIfInvalidSignature()`](https://cogdisreslab.github.io/drugfindR/reference/stopIfInvalidSignature.md)
  (indirectly used by
  [`getConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/getConcordants.md))

- Filtering parameter validation via
  [`.validateFilterSignatureInput()`](https://cogdisreslab.github.io/drugfindR/reference/dot-validateFilterSignatureInput.md)

- Concordance / network errors via internal iLINCS response processors.

If both `filterThreshold` and `filterProp` are supplied an error is
raised upstream in
[`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md).
Provide only one.

## Examples

``` r
NULL
#> NULL
```
