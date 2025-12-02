# Filter the L1000 Signature **\[stable\]**

This function filters the L1000 signature to a given threshold,
identifying up-regulated, down-regulated, or both up- and down-regulated
genes. The function supports both absolute threshold filtering and
proportional filtering based on quantiles of the expression data.

## Usage

``` r
filterSignature(signature, direction = "any", threshold = NULL, prop = NULL)
```

## Arguments

- signature:

  A data.frame, tibble, or DataFrame containing the L1000 signature.
  Must contain a column named "Value_LogDiffExp" with log fold-change
  values.

- direction:

  Character string specifying the direction to filter. Must be one of
  "up" (up-regulated genes only), "down" (down-regulated genes only), or
  "any" (both up- and down-regulated genes). Defaults to "any".

- threshold:

  Numeric value or vector specifying the log fold-change threshold(s).
  Can be: \* A single positive value: Creates symmetric thresholds
  (\\\pm threshold\\) \* A vector of two values: First value is the
  down-regulated threshold, second value is the up-regulated threshold
  Cannot be specified together with `prop`. One of `threshold` or `prop`
  must be provided.

- prop:

  Numeric value between 0 and 1 specifying the proportion of genes to
  select from the top and bottom of the expression distribution. For
  example, `prop = 0.1` selects the top 10% most up-regulated and bottom
  10% most down-regulated genes. Cannot be specified together with
  `threshold`.

## Value

A tibble containing the filtered L1000 signature with the same structure
as the input but containing only genes that meet the filtering criteria.

## Details

The filtering process follows these steps:

1.  Input validation: Checks data frame structure and parameter
    consistency

2.  Threshold calculation: Computes filtering thresholds based on either
    absolute values (`threshold`) or quantiles (`prop`)

3.  Direction-based filtering: Applies the computed thresholds according
    to the specified direction

When using `threshold`:

- Single value: Genes with \|logFC\| \>= threshold are retained

- Two values: Genes with logFC \<= `threshold[1]` OR logFC \>=
  `threshold[2]`

When using `prop`:

- Thresholds are calculated as quantiles of the expression distribution

- Down threshold = quantile(logFC, prop)

- Up threshold = quantile(logFC, 1 - prop)

## See also

`\link{getSignature}` for retrieving L1000 signatures from iLINCS,
`\link{prepareSignature}` for preparing custom signatures,
`\link{getConcordants}` for finding concordant signatures

## Examples

``` r
# Create a mock signature for demonstration
mockSignature <- data.frame(
    signatureID = rep("MOCK001", 20),
    Name_GeneSymbol = paste0("GENE", 1:20),
    ID_geneid = 1:20,
    Value_LogDiffExp = c(
        -3.5, -2.8, -2.1, -1.5, -1.2, -0.8, -0.5, -0.3,
        -0.1, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6, 2.0, 2.4, 2.9, 3.3, 3.8
    )
)

# Example 1: Filter by symmetric absolute threshold
# Keeps genes with |logFC| >= 1.5
filteredSymmetric <- filterSignature(mockSignature, threshold = 1.5)
nrow(filteredSymmetric) # Should return 8 genes
#> [1] 10

# Example 2: Filter by asymmetric absolute thresholds
# Keeps genes with logFC <= -2.0 OR logFC >= 2.5
filteredAsymmetric <- filterSignature(mockSignature, threshold = c(-2.0, 2.5))
nrow(filteredAsymmetric) # Should return 5 genes
#> [1] 6

# Example 3: Filter by proportion (top and bottom 20%)
filteredProportion <- filterSignature(mockSignature, prop = 0.2)
nrow(filteredProportion) # Should return 8 genes (4 up + 4 down)
#> [1] 8

# Example 4: Filter only up-regulated genes by threshold
upRegulated <- filterSignature(mockSignature, direction = "up", threshold = 1.0)
all(upRegulated$Value_LogDiffExp >= 1.0) # Should be TRUE
#> [1] TRUE

# Example 5: Filter only down-regulated genes by threshold
downRegulated <- filterSignature(mockSignature, direction = "down", threshold = 1.0)
all(downRegulated$Value_LogDiffExp <= -1.0) # Should be TRUE
#> [1] TRUE

# Network-dependent examples using real iLINCS data
# Get the L1000 signature for LINCSKD_28
kdSignature <- getSignature("LINCSKD_28")

# Filter for top 5% most extreme genes
topExtreme <- filterSignature(kdSignature, prop = 0.05)

# Get top 20% most up-regulated genes
topUpregulated <- filterSignature(kdSignature, direction = "up", prop = 0.2)
```
