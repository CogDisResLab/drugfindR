# Getting Started with drugfindR

## Introduction

Welcome to **drugfindR**! This package provides R-based programmatic
access to the [iLINCS](http://www.ilincs.org) (Integrative Library of
Integrated Network-Based Cellular Signatures) database for drug
repurposing and functional genomics research.

### What is drugfindR?

drugfindR enables researchers to:

- **Identify candidate repurposable drugs** from transcriptomic
  signatures
- **Discover functional relationships** between genes and compounds
- **Query LINCS databases** programmatically without web platform
  dependencies
- **Process signatures efficiently** with high-throughput batch
  capabilities

### What is LINCS?

The Library of Integrated Network-Based Cellular Signatures (LINCS)
project systematically catalogs cellular responses to genetic and
chemical perturbations. The iLINCS platform integrates three major
signature types:

- **Chemical Perturbagen (CP)**: Drug and compound treatment signatures
- **Gene Knockdown (KD)**: Gene silencing signatures
- **Gene Overexpression (OE)**: Gene upregulation signatures

## Installation

### From r-universe (Recommended)

``` r

install.packages("drugfindR",
    repos = c(
        "https://cogdisreslab.r-universe.dev",
        "https://cran.r-project.org"
    )
)
```

### From GitHub (Development Version)

``` r

if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("CogDisResLab/drugfindR")
```

## Loading the Package

``` r

library(drugfindR)
library(dplyr)
library(readr)
```

## Two Approaches to Using drugfindR

drugfindR offers two complementary approaches:

### 1. High-Level Convenience Functions

These wrapper functions handle the entire workflow with sensible
defaults:

- [`investigateSignature()`](https://cogdisreslab.github.io/drugfindR/reference/investigateSignature.md):
  Analyze transcriptomic data to find concordant drugs/genes
- [`investigateTarget()`](https://cogdisreslab.github.io/drugfindR/reference/investigateTarget.md):
  Explore functional relationships for a specific gene or drug

**Use when**: You want rapid results with minimal code

### 2. Modular Pipeline Functions

Five building-block functions provide fine-grained control:

1.  [`getSignature()`](https://cogdisreslab.github.io/drugfindR/reference/getSignature.md):
    Retrieve LINCS signatures by ID
2.  [`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.md):
    Format custom signatures for iLINCS
3.  [`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md):
    Apply thresholds to signatures
4.  [`getConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/getConcordants.md):
    Query iLINCS for concordant signatures
5.  [`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md):
    Generate consensus rankings

**Use when**: You need custom workflows or detailed intermediate results

## Quick Start Example 1: Investigate a Transcriptomic Signature

Let’s identify candidate drugs for a disease signature using the
convenience function.

### Load Example Data

We’ll use a COVID-19 differential expression dataset included with the
package:

``` r

# Load the COVID-19 differential expression data
covid_file <- system.file("extdata", "dCovid_diffexp.tsv", package = "drugfindR")
covid_diffexp <- read_tsv(covid_file, show_col_types = FALSE)

# Preview the data
head(covid_diffexp)
#> # A tibble: 6 × 3
#>   hgnc_symbol logFC     PValue
#>   <chr>       <dbl>      <dbl>
#> 1 CCL4L2      -3.98 0.00000177
#> 2 IL5RA       -4.83 0.00000870
#> 3 FN1          4.94 0.0000117 
#> 4 GSTM1       -8.21 0.0000153 
#> 5 CD180        2.15 0.0000202 
#> 6 FAM20C       3.11 0.0000255
```

This dataset contains: - `hgnc_symbol`: Gene symbols - `logFC`: Log
fold-change values - `PValue`: Statistical significance

### One-Line Analysis

``` r

# Find drugs that may counteract the COVID-19 signature
results <- investigateSignature(
    expr = covid_diffexp,
    outputLib = "CP", # Query Chemical Perturbagens
    filterThreshold = 1.5, # Keep genes with |logFC| >= 1.5
    geneColumn = "hgnc_symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)

# View top candidates
head(results, 10)
```

### Understanding the Results

The output contains:

- **Source**: Your input signature identifier
- **Target**: Drug/compound name
- **Similarity**: Concordance score (-1 to 1)
  - Negative scores: Drug reverses your signature (potential
    therapeutic)
  - Positive scores: Drug mimics your signature
- **TargetSignature**: iLINCS signature ID
- **TargetCellLine**: Cell line used in the signature
- **pValue**: Statistical significance

## Quick Start Example 2: Investigate a Specific Gene

Let’s explore what happens when TP53 is knocked down and find drugs with
similar effects.

``` r

# Find drugs that mimic TP53 knockdown
tp53_results <- investigateTarget(
    target = "TP53",
    inputLib = "KD", # Use Knockdown signatures
    outputLib = "CP", # Find Chemical Perturbagens
    filterThreshold = 0.5
)

# View results
head(tp53_results, 10)
```

## Modular Approach: Step-by-Step Workflow

For more control, let’s break down the analysis into individual steps.

### Step 1: Prepare Your Signature

``` r

# Convert your data to iLINCS format
signature <- prepareSignature(
    covid_diffexp,
    geneColumn = "hgnc_symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)

# Preview the prepared signature
head(signature)
#>    signatureID ID_geneid Name_GeneSymbol Value_LogDiffExp Significance_pvalue
#> 1     InputSig      4860             PNP         1.709692         0.002436390
#> 5     InputSig      4125          MAN2B1         1.100506         0.003027542
#> 8     InputSig      2274            FHL2         1.330287         0.142375657
#> 14    InputSig       351             APP         1.050513         0.010844053
#> 26    InputSig      7077           TIMP2         1.795990         0.012416655
#> 29    InputSig     23659         PLA2G15         1.113302         0.049536693
```

The
[`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.md)
function: - Standardizes column names - Maps genes to L1000 gene space -
Validates data format

### Step 2: Filter by Threshold

``` r

# Filter for strongly upregulated genes
filtered_up <- filterSignature(
    signature,
    direction = "up",
    threshold = 1.5
)

# Filter for strongly downregulated genes
filtered_down <- filterSignature(
    signature,
    direction = "down",
    threshold = 1.5
)

cat("Upregulated genes:", nrow(filtered_up), "\n")
#> Upregulated genes: 72
cat("Downregulated genes:", nrow(filtered_down), "\n")
#> Downregulated genes: 15
```

### Step 3: Query for Concordant Signatures

``` r

# Find drugs that match upregulated genes
concordants_up <- getConcordants(
    filtered_up,
    ilincsLibrary = "CP",
    direction = "up"
)

# Find drugs that match downregulated genes
concordants_down <- getConcordants(
    filtered_down,
    ilincsLibrary = "CP",
    direction = "down"
)
```

### Step 4: Generate Consensus Rankings

``` r

# Combine and rank candidate drugs
consensus <- consensusConcordants(
    concordants_up,
    concordants_down,
    paired = TRUE,
    cutoff = 0.2
)

head(consensus, 10)
```

## Choosing Your Approach

| Aspect | Convenience Functions | Modular Functions |
|----|----|----|
| **Code Length** | Minimal (1 function call) | Longer (5 steps) |
| **Flexibility** | Limited customization | Full control |
| **Intermediate Results** | Hidden | Available |
| **Learning Curve** | Easy | Moderate |
| **Best For** | Quick analyses, standard workflows | Custom pipelines, research |

## Filtering Strategies

### Absolute Thresholds

Use specific log fold-change cutoffs:

``` r

# Single threshold (symmetric)
filtered_symmetric <- filterSignature(signature, threshold = 1.5)
# Keeps genes with |logFC| >= 1.5

# Double threshold (asymmetric)
filtered_asymmetric <- filterSignature(signature, threshold = c(-2.0, 1.5))
# Keeps genes with logFC <= -2.0 OR logFC >= 1.5
```

### Proportional Thresholds

Select a percentage of most extreme genes:

``` r

# Top/bottom 10% most differentially expressed
filtered_prop <- filterSignature(signature, prop = 0.1)
```

### Direction-Specific Filtering

``` r

# Only upregulated genes
up_only <- filterSignature(signature, direction = "up", threshold = 1.0)

# Only downregulated genes
down_only <- filterSignature(signature, direction = "down", threshold = 1.0)

# Both directions (default)
both_directions <- filterSignature(signature, direction = "any", threshold = 1.0)
```

## Understanding Library Types

### Chemical Perturbagen (CP)

- Drug and compound treatment signatures
- Use for: Drug repurposing, finding therapeutic candidates
- Example: `outputLib = "CP"`

### Gene Knockdown (KD)

- Gene silencing (shRNA, siRNA) signatures
- Use for: Understanding gene function, finding genes with similar
  effects
- Example: `inputLib = "KD"`

### Gene Overexpression (OE)

- Gene upregulation signatures
- Use for: Understanding gain-of-function effects
- Example: `inputLib = "OE"`

## Paired vs. Unpaired Analysis

### Paired Analysis (Default)

Separately analyzes up- and down-regulated genes:

``` r

# Default behavior
results_paired <- investigateSignature(
    covid_diffexp,
    outputLib = "CP",
    filterThreshold = 1.5,
    paired = TRUE, # This is the default
    geneColumn = "hgnc_symbol"
)
```

**Advantages**: - More biological precision - Distinguishes directional
effects - Better for complex signatures

### Unpaired Analysis

Treats all significant genes together:

``` r

results_unpaired <- investigateSignature(
    covid_diffexp,
    outputLib = "CP",
    filterThreshold = 1.5,
    paired = FALSE,
    geneColumn = "hgnc_symbol"
)
```

**Advantages**: - Simpler interpretation - Faster execution - Good for
exploratory analysis

## Common Use Cases

### 1. Drug Repurposing

**Goal**: Find drugs that reverse a disease signature

``` r

disease_results <- investigateSignature(
    expr = disease_signature,
    outputLib = "CP",
    filterThreshold = 1.5,
    paired = TRUE
)

# Look for negative similarity scores (drug reverses disease)
potential_therapeutics <- disease_results %>%
    filter(Similarity < -0.3) %>%
    arrange(Similarity)
```

### 2. Gene Function Discovery

**Goal**: Understand what a gene does by finding similar genes

``` r

gene_function <- investigateTarget(
    target = "BRCA1",
    inputLib = "KD",
    outputLib = "KD", # Find genes with similar knockdown effects
    filterThreshold = 0.5
)
```

### 3. Mechanism of Action

**Goal**: Find genes affected by a drug

``` r

drug_mechanism <- investigateTarget(
    target = "metformin",
    inputLib = "CP",
    outputLib = "KD", # Find genes whose KD mimics drug effect
    filterThreshold = 0.5
)
```

## Next Steps

Now that you understand the basics, explore these advanced topics:

- **[Drug Repurposing from Transcriptomic
  Data](https://cogdisreslab.github.io/drugfindR/articles/drug-repurposing.md)**:
  Deep dive into therapeutic discovery workflows
- **[Target Investigation and Functional
  Genomics](https://cogdisreslab.github.io/drugfindR/articles/target-investigation.md)**:
  Comprehensive guide to gene and drug target analysis
- **[Package
  Documentation](https://cogdisreslab.github.io/drugfindR/)**: Complete
  function reference

## Session Information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] readr_2.2.0         dplyr_1.2.1         drugfindR_0.99.1170
#> [4] BiocStyle_2.40.0   
#> 
#> loaded via a namespace (and not attached):
#>  [1] utf8_1.2.6          rappdirs_0.3.4      sass_0.4.10        
#>  [4] generics_0.1.4      stringi_1.8.7       hms_1.1.4          
#>  [7] digest_0.6.39       magrittr_2.0.5      evaluate_1.0.5     
#> [10] bookdown_0.47       fastmap_1.2.0       jsonlite_2.0.0     
#> [13] DFplyr_1.6.0        BiocManager_1.30.27 purrr_1.2.2        
#> [16] httr2_1.2.3         textshaping_1.0.5   jquerylib_0.1.4    
#> [19] cli_3.6.6           crayon_1.5.3        rlang_1.3.0        
#> [22] bit64_4.8.2         withr_3.0.3         cachem_1.1.0       
#> [25] yaml_2.3.12         otel_0.2.0          parallel_4.6.1     
#> [28] tools_4.6.1         tzdb_0.5.0          BiocGenerics_0.58.1
#> [31] curl_7.1.0          vctrs_0.7.3         R6_2.6.1           
#> [34] stats4_4.6.1        lifecycle_1.0.5     stringr_1.6.0      
#> [37] S4Vectors_0.50.1    fs_2.1.0            htmlwidgets_1.6.4  
#> [40] bit_4.6.0           vroom_1.7.1         ragg_1.5.2         
#> [43] pkgconfig_2.0.3     desc_1.4.3          pkgdown_2.2.1      
#> [46] pillar_1.11.1       bslib_0.11.0        glue_1.8.1         
#> [49] systemfonts_1.3.2   xfun_0.60           tibble_3.3.1       
#> [52] tidyselect_1.2.1    knitr_1.51          htmltools_0.5.9    
#> [55] rmarkdown_2.31      compiler_4.6.1
```

## References

1.  O’Donovan SM, Imami A, et al. (2021). Identification of candidate
    repurposable drugs to combat COVID-19 using a signature-based
    approach. *Scientific Reports*, 11:4495.
    [doi:10.1038/s41598-021-84044-9](https://doi.org/10.1038/s41598-021-84044-9)

2.  Subramanian A, et al. (2017). A Next Generation Connectivity Map:
    L1000 Platform and the First 1,000,000 Profiles. *Cell*,
    171(6):1437-1452.e17.
    [doi:10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)
