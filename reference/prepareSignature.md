# Prepare an L1000 Signature from a given differential gene expression output **\[stable\]**

This function takes a differential gene expression output from any
pipeline like edgeR or DeSeq2 or any that give you the gene symbol,
log_2 fold-change and p-value and transforms that into an L1000
signature for later processing.

## Usage

``` r
prepareSignature(
  dge,
  geneColumn = "Symbol",
  logfcColumn = "logFC",
  pvalColumn = "PValue"
)
```

## Arguments

- dge:

  A dataframe-like object that has the differential gene expression
  information

- geneColumn:

  The name of the column that has gene symbols

- logfcColumn:

  The name of the column that has log_2 fold-change values

- pvalColumn:

  The name of the column that has p-values

## Value

A tibble with the L1000 signature.

## Examples

``` r
# Load example differential expression data from package
dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
    package = "drugfindR"
)
dge_data <- read.delim(dge_file)

# Prepare signature with p-values (standard workflow)
signature <- prepareSignature(
    dge_data,
    geneColumn = "hgnc_symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)
head(signature)
#>    signatureID ID_geneid Name_GeneSymbol Value_LogDiffExp Significance_pvalue
#> 1     InputSig      4860             PNP         1.709692         0.002436390
#> 5     InputSig      4125          MAN2B1         1.100506         0.003027542
#> 8     InputSig      2274            FHL2         1.330287         0.142375657
#> 14    InputSig       351             APP         1.050513         0.010844053
#> 26    InputSig      7077           TIMP2         1.795990         0.012416655
#> 29    InputSig     23659         PLA2G15         1.113302         0.049536693

# Prepare signature without p-values
signature_no_pval <- prepareSignature(
    dge_data,
    geneColumn = "hgnc_symbol",
    logfcColumn = "logFC",
    pvalColumn = NA
)
head(signature_no_pval)
#>    signatureID ID_geneid Name_GeneSymbol Value_LogDiffExp
#> 1     InputSig      4860             PNP         1.709692
#> 5     InputSig      4125          MAN2B1         1.100506
#> 8     InputSig      2274            FHL2         1.330287
#> 14    InputSig       351             APP         1.050513
#> 26    InputSig      7077           TIMP2         1.795990
#> 29    InputSig     23659         PLA2G15         1.113302

# Custom column names example
custom_dge <- data.frame(
    Gene = c("TP53", "MYC", "BRCA1", "EGFR"),
    FC = c(2.5, -1.8, 3.2, -2.1),
    Pval = c(0.001, 0.01, 0.0001, 0.005)
)

custom_signature <- prepareSignature(
    custom_dge,
    geneColumn = "Gene",
    logfcColumn = "FC",
    pvalColumn = "Pval"
)
print(custom_signature)
#>    signatureID ID_geneid Name_GeneSymbol Value_LogDiffExp Significance_pvalue
#> 1     InputSig      1956            EGFR             -2.1               5e-03
#> 10    InputSig      7157            TP53              2.5               1e-03
#> 16    InputSig       672           BRCA1              3.2               1e-04
#> 26    InputSig      4609             MYC             -1.8               1e-02
```
