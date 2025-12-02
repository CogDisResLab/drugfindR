# Rename target-related columns to user-facing output names

Standardizes internal concordants/consensus column names to the
user-facing output schema used by this package.

## Usage

``` r
targetRename(inputNames)
```

## Arguments

- inputNames:

  Character vector of column names to rename. See Details for the
  expected input ordering and mapping.

## Value

Character vector of output (renamed) column names.

## Details

Expected input columns (by position) are the internal concordants
fields:

1.  `signatureid`, 2) `treatment` (or `compound` pre-renaming),

2.  `cellline`, 4) `time`, 5) `concentration`, 6) `sig_direction`,

3.  `sig_type`, 8) `similarity`, 9) `pValue`.

These are mapped to the user-facing names returned by functions like
[`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md)
and downstream investigation helpers:

- `TargetSignature`, `Target`, `TargetCellLine`, `TargetTime`,
  `TargetConcentration`, `InputSigDirection`, `SignatureType`,
  `Similarity`, `pValue`.

Only the names are returned; the renaming is applied via
`dplyr::rename_with(x, targetRename)`.

## Examples

``` r
NULL
#> NULL
```
