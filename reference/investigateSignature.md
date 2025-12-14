# Investigate a given DGE dataset **\[stable\]**

This function takes a DGE Data frame and then finds concordant
signatures to that. This generates an L1000 signature from the DGE
dataset and then uploads that signature to iLINCS to find the relevant
concordant (or discordant) signatures

## Usage

``` r
investigateSignature(
  expr,
  outputLib,
  filterThreshold = NULL,
  filterProp = NULL,
  similarityThreshold = 0.2,
  paired = TRUE,
  outputCellLines = NULL,
  geneColumn = "Symbol",
  logfcColumn = "logFC",
  pvalColumn = "PValue",
  sourceName = "Input",
  sourceCellLine = NA,
  sourceTime = NA,
  sourceConcentration = NA
)
```

## Arguments

- expr:

  A dataframe that has differential gene expression analysis

- outputLib:

  The library to search

- filterThreshold:

  The Filtering threshold.

- filterProp:

  The Filtering proportion.

- similarityThreshold:

  The Similarity Threshold

- paired:

  Logical. Whether to query iLINCS separately for up and down regulated
  genes

- outputCellLines:

  A character vector of cell lines to restrict the output search to.

- geneColumn:

  The name of the column that has gene symbols

- logfcColumn:

  The name of the column that has log_2 fold-change values

- pvalColumn:

  The name of the column that has p-values

- sourceName:

  (Optional) An annotation column to identify the signature by name

- sourceCellLine:

  (Optional) An annotation column to specify the cell line for the input
  data

- sourceTime:

  (Optional) An annotation column to specify the time for the input data

- sourceConcentration:

  (Optional) An annotation column to specify the concentration for the
  input data

## Value

A tibble with the the similarity scores and signature metadata

## Examples

``` r
# Input validation example (no API calls)
mockExpr <- data.frame(
    Symbol = c("TP53", "MYC"),
    logFC = c(2.5, -1.8),
    PValue = c(0.001, 0.01)
)

# Validate library parameter (should produce error)
tryCatch(
    investigateSignature(mockExpr, outputLib = "INVALID"),
    error = function(e) message("Expected error: invalid library")
)
#> Expected error: invalid library

# \donttest{
# This function makes multiple API calls to iLINCS and may take several minutes

# Load differential expression data
inputSignature <- read.table(
    system.file("extdata", "dCovid_diffexp.tsv", package = "drugfindR"),
    header = TRUE
)

# Investigate the signature against chemical perturbagen library
investigatedSignature <- investigateSignature(
    inputSignature,
    outputLib = "CP",
    filterThreshold = 0.5,
    geneColumn = "hgnc_symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)
head(investigatedSignature)
#> # A tibble: 6 × 14
#>   Source Target    Similarity SourceSignature SourceCellLine SourceConcentration
#>   <chr>  <chr>          <dbl> <chr>           <lgl>          <lgl>              
#> 1 Input  FR-180204     -0.760 InputSig        NA             NA                 
#> 2 Input  Ponatinib      0.751 InputSig        NA             NA                 
#> 3 Input  Suprofen       0.740 InputSig        NA             NA                 
#> 4 Input  SGI-1776       0.739 InputSig        NA             NA                 
#> 5 Input  SCHEMBL6…     -0.733 InputSig        NA             NA                 
#> 6 Input  Lypressin     -0.729 InputSig        NA             NA                 
#> # ℹ 8 more variables: SourceTime <lgl>, TargetSignature <chr>,
#> #   TargetCellLine <chr>, TargetConcentration <chr>, TargetTime <chr>,
#> #   InputSigDirection <chr>, SignatureType <chr>, pValue <dbl>
# }
```
