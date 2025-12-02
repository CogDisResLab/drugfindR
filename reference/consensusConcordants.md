# Generate a Consensus list of Targets **\[stable\]**

This function takes a list of (optionally split) concordance dataframes
and returns a ranked list of gene or drug targets that have been chose
for their maximal similarity to the signature

## Usage

``` r
consensusConcordants(..., paired = FALSE, cutoff = 0.321, cellLine = NULL)
```

## Arguments

- ...:

  One or Two (see paired) Data Frames with the concordants

- paired:

  Logical indicating whether you split the dataframes by up and down
  regulated in prior analysis

- cutoff:

  A similarity cutoff value. Defaults to 0.321

- cellLine:

  A character vector of Cell Lines you are interested in.

## Value

A tibble with the filtered and deduplicated results

## Examples

``` r
# Create mock concordants data for demonstration
mockConcordants <- data.frame(
    signatureid = paste0("SIG", 1:10),
    treatment = c(
        "TP53", "TP53", "MYC", "MYC", "EGFR",
        "EGFR", "KRAS", "BRCA1", "BRCA1", "PIK3CA"
    ),
    cellline = c(
        "A375", "PC3", "A375", "MCF7", "A375",
        "PC3", "A375", "A375", "MCF7", "A375"
    ),
    time = rep("24H", 10),
    concentration = rep(NA, 10),
    sig_direction = rep("DOWN", 10),
    sig_type = rep("single", 10),
    similarity = c(
        0.85, 0.72, -0.68, -0.45, 0.55,
        0.38, 0.42, 0.51, 0.33, 0.29
    ),
    pValue = rep(0.001, 10)
)

# Example 1: Basic consensus with default cutoff
consensus <- consensusConcordants(mockConcordants)
nrow(consensus) # Targets with |similarity| >= 0.321
#> [1] 5

# Example 2: Consensus with higher cutoff
consensus_strict <- consensusConcordants(mockConcordants, cutoff = 0.5)
nrow(consensus_strict) # Fewer targets with higher threshold
#> [1] 4

# Example 3: Filter by cell line
consensus_A375 <- consensusConcordants(mockConcordants, cellLine = "A375")
unique(consensus_A375$CellLine) # Only A375
#> Warning: Unknown or uninitialised column: `CellLine`.
#> NULL

# \donttest{
# Network-dependent examples using real iLINCS data
# Get the L1000 signature for LINCSKD_28
kdSignature <- getSignature("LINCSKD_28")

# Get concordant gene knockdown signatures
concordantSignatures <- getConcordants(kdSignature, ilincsLibrary = "KD")

# Get the consensus list with different parameters
consensus <- consensusConcordants(concordantSignatures, cutoff = 0.5)

# Paired analysis example
filteredUp <- filterSignature(kdSignature, direction = "up", threshold = 0.5)
filteredDown <- filterSignature(kdSignature, direction = "down", threshold = -0.5)
concordants_up <- getConcordants(filteredUp, ilincsLibrary = "KD")
concordants_down <- getConcordants(filteredDown, ilincsLibrary = "KD")
consensus <- consensusConcordants(concordants_up, concordants_down, paired = TRUE)
# }
```
