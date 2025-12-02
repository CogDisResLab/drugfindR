# Apply filtering based on direction and thresholds

This internal function performs the actual filtering of the signature
data based on the specified direction and calculated thresholds. It
implements the core filtering logic using dplyr operations.

## Usage

``` r
.applyDirectionFilter(signature, thresholds, direction = "any")
```

## Arguments

- signature:

  A data.frame-like object containing the signature data. Must have a
  column named "Value_LogDiffExp" containing log fold-change values.

- thresholds:

  A named list containing:

  - `downThreshold`: Threshold for down-regulated genes

  - `upThreshold`: Threshold for up-regulated genes

- direction:

  Character string specifying the filtering direction. Must be one of:
  \* "up": Keep only up-regulated genes (logFC \>= upThreshold) \*
  "down": Keep only down-regulated genes (logFC \<= downThreshold) \*
  "any": Keep both up- and down-regulated genes (logFC \>= upThreshold
  OR logFC \<= downThreshold)

## Value

A tibble containing the filtered signature data with the same structure
as the input but including only rows that meet the filtering criteria.

## Details

The filtering logic depends on the direction parameter:

- `"up"`: Retains genes where `Value_LogDiffExp >= upThreshold`

- `"down"`: Retains genes where `Value_LogDiffExp <= downThreshold`

- `"any"`: Retains genes where
  `Value_LogDiffExp >= upThreshold OR Value_LogDiffExp <= downThreshold`

The function uses
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
with [`rlang::.data`](https://rlang.r-lib.org/reference/dot-data.html)
for non-standard evaluation, ensuring compatibility with different data
frame types and avoiding issues with variable scoping.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create sample signature data
signature <- data.frame(
    signatureID = rep("TEST", 10),
    Name_GeneSymbol = paste0("GENE", 1:10),
    Value_LogDiffExp = c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4)
)

# Define thresholds
thresholds <- list(downThreshold = -1.5, upThreshold = 1.5)

# Filter for up-regulated genes only
up_filtered <- .applyDirectionFilter(signature, "up", thresholds)
# Returns genes with logFC >= 1.5 (GENE8, GENE9, GENE10)

# Filter for down-regulated genes only
down_filtered <- .applyDirectionFilter(signature, "down", thresholds)
# Returns genes with logFC <= -1.5 (GENE1, GENE2)

# Filter for both up- and down-regulated genes
both_filtered <- .applyDirectionFilter(signature, "any", thresholds)
# Returns genes with |logFC| >= 1.5 (GENE1, GENE2, GENE8, GENE9, GENE10)
} # }
```
