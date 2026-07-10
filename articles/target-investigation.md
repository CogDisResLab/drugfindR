# Target Investigation and Functional Genomics

## Introduction

The
[`investigateTarget()`](https://cogdisreslab.github.io/drugfindR/reference/investigateTarget.md)
function in **drugfindR** enables comprehensive functional genomics
analyses by exploring relationships between genes and drugs in the LINCS
database. This vignette covers:

- Understanding what happens when a gene is knocked down or
  overexpressed
- Finding drugs that mimic or reverse gene perturbation effects
- Discovering genes with similar functional roles
- Elucidating drug mechanisms of action

## Prerequisites

``` r

library(drugfindR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
```

## The Target Investigation Workflow

### Key Questions Answered

1.  **Gene Function**: What cellular changes occur when a gene is
    perturbed?
2.  **Drug Mimicry**: Which drugs produce effects similar to gene
    perturbation?
3.  **Drug Opposition**: Which drugs counteract gene perturbation
    effects?
4.  **Functional Similarity**: Which genes have similar functional
    roles?
5.  **Mechanism of Action**: What genes does a drug affect?

### Workflow Overview

    Target (Gene/Drug)
        ↓
    Retrieve LINCS Signatures
        ↓
    Filter by Magnitude/Direction
        ↓
    Query Output Library
        ↓
    Generate Consensus
        ↓
    Functional Insights

## Basic Target Investigation

### Example 1: What does TP53 knockdown do?

TP53 is a tumor suppressor gene. Let’s explore what happens when it’s
knocked down:

``` r

# Investigate TP53 knockdown effects
tp53_kd <- investigateTarget(
    target = "TP53",
    inputLib = "KD", # Use knockdown signatures
    outputLib = "KD", # Find similar gene knockdowns
    filterThreshold = 0.5,
    similarityThreshold = 0.3
)

# View results
head(tp53_kd, 15)
```

#### Interpreting Results

The output shows:

- **Target**: Other genes whose knockdown produces similar effects
- **Similarity**: Positive values = similar function, negative =
  opposite function
- **SourceSignature**: TP53 knockdown signature ID(s)
- **TargetSignature**: Matched gene signature IDs

``` r

# Find genes with similar function to TP53
similar_genes <- tp53_kd %>%
    filter(Similarity > 0.3) %>%
    arrange(desc(Similarity))

cat("Genes functionally similar to TP53:\n")
head(similar_genes, 10)

# Find genes with opposite function
opposite_genes <- tp53_kd %>%
    filter(Similarity < -0.3) %>%
    arrange(Similarity)

cat("\nGenes with opposite function:\n")
head(opposite_genes, 10)
```

### Example 2: Which drugs mimic TP53 loss?

Cancer cells often lose TP53 function. Which drugs produce similar
effects?

``` r

# Find drugs that mimic TP53 knockdown
tp53_mimics <- investigateTarget(
    target = "TP53",
    inputLib = "KD", # Start with KD signature
    outputLib = "CP", # Find chemical perturbagens
    filterThreshold = 0.5,
    similarityThreshold = 0.3
)

# Drugs with positive similarity mimic TP53 loss
mimicking_drugs <- tp53_mimics %>%
    filter(Similarity > 0.3) %>%
    arrange(desc(Similarity))

head(mimicking_drugs, 10)
```

### Example 3: Which drugs rescue TP53 loss?

More importantly, which drugs might compensate for TP53 loss?

``` r

# Drugs with negative similarity oppose TP53 loss
rescue_drugs <- tp53_mimics %>%
    filter(Similarity < -0.3) %>%
    arrange(Similarity)

cat("Potential TP53 loss rescue compounds:\n")
head(rescue_drugs, 15)
```

## Gene Overexpression Analysis

### Example 4: Effects of MYC overexpression

MYC is an oncogene frequently overexpressed in cancer:

``` r

# Investigate MYC overexpression
myc_oe <- investigateTarget(
    target = "MYC",
    inputLib = "OE", # Overexpression signatures
    outputLib = "CP", # Find drugs
    filterThreshold = 0.5
)

# Drugs that oppose MYC overexpression (therapeutic potential)
myc_antagonists <- myc_oe %>%
    filter(Similarity < -0.3) %>%
    arrange(Similarity)

cat("Drugs that oppose MYC overexpression:\n")
head(myc_antagonists, 10)
```

### Comparing KD vs OE

Different perturbation types reveal different biology:

``` r

# BRCA1 knockdown
brca1_kd <- investigateTarget(
    target = "BRCA1",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5
)

# BRCA1 overexpression
brca1_oe <- investigateTarget(
    target = "BRCA1",
    inputLib = "OE",
    outputLib = "CP",
    filterThreshold = 0.5
)

# Compare top drugs
kd_top <- brca1_kd %>%
    filter(Similarity < -0.3) %>%
    pull(Target)
oe_top <- brca1_oe %>%
    filter(Similarity < -0.3) %>%
    pull(Target)

overlap <- intersect(kd_top, oe_top)
cat("Drugs opposing both KD and OE:", length(overlap), "\n")
```

## Drug Mechanism of Action

### Example 5: What does metformin affect?

Reverse the query to understand drug mechanisms:

``` r

# Find genes affected by metformin
metformin_targets <- investigateTarget(
    target = "metformin",
    inputLib = "CP", # Drug signature
    outputLib = "KD", # Find similar gene knockdowns
    filterThreshold = 0.5
)

# Genes whose loss mimics metformin treatment
affected_genes <- metformin_targets %>%
    filter(Similarity > 0.3) %>%
    arrange(desc(Similarity))

cat("Genes potentially targeted by metformin:\n")
head(affected_genes, 15)
```

### Example 6: Compare drug to drug

Find drugs with similar mechanisms:

``` r

# Find drugs similar to aspirin
aspirin_similar <- investigateTarget(
    target = "aspirin",
    inputLib = "CP",
    outputLib = "CP",
    filterThreshold = 0.5
)

similar_drugs <- aspirin_similar %>%
    filter(Similarity > 0.4) %>%
    arrange(desc(Similarity))

head(similar_drugs, 10)
```

## Paired vs. Unpaired Analysis

### Paired Analysis (Default)

Analyzes up and down-regulated genes separately:

``` r

# Paired analysis captures directional effects
egfr_paired <- investigateTarget(
    target = "EGFR",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5,
    paired = TRUE # Default
)

# Results preserve directional information
cat("Directional analysis:\n")
table(egfr_paired$InputSigDirection)
```

#### When to use paired:

- Complex genes with many targets
- When directionality matters
- More biologically precise
- Recommended for most analyses

### Unpaired Analysis

Combines all significant genes:

``` r

# Unpaired analysis for simpler interpretation
egfr_unpaired <- investigateTarget(
    target = "EGFR",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5,
    paired = FALSE
)

# Faster, simpler results
cat("Combined analysis genes:\n")
nrow(egfr_unpaired)
```

#### When to use unpaired:

- Exploratory analysis
- Simpler interpretation needed
- Symmetric signatures
- Computational efficiency matters

## Cell Line Considerations

### Filtering Input Cell Lines

Restrict source signatures to specific cell lines:

``` r

# Only use TP53 KD from breast cancer cell lines
tp53_breast <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    inputCellLines = c("MCF7", "MDAMB231", "T47D"),
    filterThreshold = 0.5
)

cat("Results from breast cancer cell lines only\n")
```

### Filtering Output Cell Lines

Restrict target signatures to specific contexts:

``` r

# Find drugs tested in lung cancer cell lines
kras_lung <- investigateTarget(
    target = "KRAS",
    inputLib = "KD",
    outputLib = "CP",
    outputCellLines = c("A549", "H1299", "H460"),
    filterThreshold = 0.5
)

cat("Drugs tested in lung cancer models\n")
```

### Cross-Cell Line Analysis

``` r

# Get all results
all_cells <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5
)

# Analyze cell line distribution
cellline_summary <- all_cells %>%
    group_by(Target) %>%
    summarize(
        n_celllines = n_distinct(TargetCellLine),
        mean_similarity = mean(Similarity),
        .groups = "drop"
    ) %>%
    filter(n_celllines >= 3) %>%
    arrange(desc(n_celllines))

cat("Drugs tested across multiple cell lines:\n")
head(cellline_summary, 10)
```

## Threshold Optimization

### Conservative Analysis

High confidence, fewer hits:

``` r

conservative <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.85, # Very stringent
    similarityThreshold = 0.4 # High similarity required
)

cat("Conservative results:", nrow(conservative), "\n")
```

### Liberal Analysis

Broader coverage, more hypotheses:

``` r

liberal <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.3, # Permissive
    similarityThreshold = 0.15 # Lower cutoff
)

cat("Liberal results:", nrow(liberal), "\n")
```

### Comparing Thresholds

``` r

# Test multiple thresholds
thresholds <- c(0.3, 0.5, 0.7, 0.85)

threshold_results <- map_dfr(thresholds, function(thresh) {
    results <- investigateTarget(
        target = "TP53",
        inputLib = "KD",
        outputLib = "CP",
        filterThreshold = thresh,
        similarityThreshold = 0.2
    )

    tibble(
        threshold = thresh,
        n_results = nrow(results),
        n_strong = sum(abs(results$Similarity) > 0.4)
    )
})

print(threshold_results)
```

## Multiple Target Analysis

### Batch Processing Targets

``` r

# Analyze multiple related genes
cell_cycle_genes <- c("TP53", "RB1", "CDKN2A", "MYC", "E2F1")

cell_cycle_results <- map_dfr(cell_cycle_genes, function(gene) {
    tryCatch(
        {
            results <- investigateTarget(
                target = gene,
                inputLib = "KD",
                outputLib = "CP",
                filterThreshold = 0.5
            )

            results %>%
                mutate(SourceGene = gene) %>%
                filter(Similarity < -0.3) # Therapeutic candidates
        },
        error = function(e) {
            message("No results for ", gene)
            tibble()
        }
    )
})

# Find drugs targeting multiple cell cycle genes
multi_target_drugs <- cell_cycle_results %>%
    group_by(Target) %>%
    summarize(
        n_genes_affected = n_distinct(SourceGene),
        mean_similarity = mean(Similarity),
        .groups = "drop"
    ) %>%
    filter(n_genes_affected >= 3) %>%
    arrange(n_genes_affected, mean_similarity)

cat("Drugs affecting multiple cell cycle genes:\n")
head(multi_target_drugs, 10)
```

## Pathway-Level Analysis

### Analyzing Gene Sets

``` r

# DNA repair pathway genes
dna_repair_genes <- c("BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2")

# Get concordant drugs for each gene
dna_repair_drugs <- map_dfr(dna_repair_genes, function(gene) {
    tryCatch(
        {
            investigateTarget(
                target = gene,
                inputLib = "KD",
                outputLib = "CP",
                filterThreshold = 0.5
            ) %>%
                mutate(PathwayGene = gene)
        },
        error = function(e) tibble()
    )
})

# Drugs consistently affecting the pathway
pathway_drugs <- dna_repair_drugs %>%
    group_by(Target) %>%
    summarize(
        genes_affected = n_distinct(PathwayGene),
        mean_effect = mean(Similarity),
        .groups = "drop"
    ) %>%
    filter(genes_affected >= 4) %>%
    arrange(desc(genes_affected))

cat("Drugs affecting DNA repair pathway:\n")
head(pathway_drugs, 10)
```

## Visualization Strategies

### Similarity Score Distribution

``` r

tp53_results <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5
)

ggplot(tp53_results, aes(x = Similarity)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(
        xintercept = c(-0.3, 0.3),
        linetype = "dashed",
        color = "red",
        size = 1
    ) +
    labs(
        title = "TP53 Knockdown: Drug Similarity Distribution",
        x = "Similarity Score",
        y = "Number of Drug Signatures",
        subtitle = "Negative = opposes TP53 loss; Positive = mimics TP53 loss"
    ) +
    theme_minimal()
```

### Top Candidates Bar Plot

``` r

top_drugs <- tp53_results %>%
    arrange(Similarity) %>%
    head(20)

ggplot(top_drugs, aes(x = reorder(Target, Similarity), y = Similarity)) +
    geom_col(aes(fill = Similarity < 0), show.legend = FALSE) +
    scale_fill_manual(values = c("red3", "steelblue")) +
    coord_flip() +
    labs(
        title = "Top 20 Drugs Related to TP53 Knockdown",
        subtitle = "Blue = opposes TP53 loss (therapeutic potential)",
        x = "Drug/Compound",
        y = "Similarity Score"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
```

### Network Visualization Concept

``` r

# Prepare data for network visualization
# (requires igraph or visNetwork package)

# Get relationships between multiple genes and drugs
network_data <- map_dfr(c("TP53", "BRCA1", "EGFR"), function(gene) {
    investigateTarget(
        target = gene,
        inputLib = "KD",
        outputLib = "CP",
        filterThreshold = 0.5,
        similarityThreshold = 0.3
    ) %>%
        mutate(GeneSource = gene) %>%
        select(GeneSource, Target, Similarity) %>%
        filter(abs(Similarity) > 0.4)
})

# This data can be used with igraph to create gene-drug networks
head(network_data)
```

### Heatmap of Gene-Drug Relationships

``` r

# Create a matrix of gene-drug similarities
genes_of_interest <- c("TP53", "BRCA1", "MYC", "EGFR", "KRAS")

heatmap_data <- map_dfr(genes_of_interest, function(gene) {
    investigateTarget(
        target = gene,
        inputLib = "KD",
        outputLib = "CP",
        filterThreshold = 0.5
    ) %>%
        mutate(Gene = gene) %>%
        group_by(Gene, Target) %>%
        summarize(Similarity = mean(Similarity), .groups = "drop")
})

# Filter to drugs affecting multiple genes
common_drugs <- heatmap_data %>%
    group_by(Target) %>%
    filter(n() >= 3) %>%
    ungroup()

ggplot(common_drugs, aes(x = Target, y = Gene, fill = Similarity)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0,
        name = "Similarity"
    ) +
    labs(
        title = "Gene-Drug Relationship Heatmap",
        x = "Drug",
        y = "Gene"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
```

## Advanced Use Cases

### Case 1: Synthetic Lethality Discovery

Find drugs that are specifically toxic when a gene is lost:

``` r

# In BRCA1-deficient cells, PARP inhibitors are synthetically lethal
# Find similar relationships

brca1_synthetic <- investigateTarget(
    target = "BRCA1",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.6,
    paired = TRUE
)

# Drugs with strong positive similarity might indicate synergy
potential_synthetic <- brca1_synthetic %>%
    filter(Similarity > 0.4) %>%
    arrange(desc(Similarity))

cat("Potential synthetic lethal partners:\n")
head(potential_synthetic, 10)
```

### Case 2: Rescue vs. Enhancement

``` r

results <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5
)

# Classify drugs
drug_classification <- results %>%
    mutate(
        Classification = case_when(
            Similarity < -0.4 ~ "Strong Rescue",
            Similarity < -0.2 ~ "Moderate Rescue",
            Similarity > 0.4 ~ "Strong Enhancement",
            Similarity > 0.2 ~ "Moderate Enhancement",
            TRUE ~ "Neutral"
        )
    )

table(drug_classification$Classification)
```

### Case 3: Temporal Analysis

``` r

# Analyze time-dependent effects
results <- investigateTarget(
    target = "MYC",
    inputLib = "OE",
    outputLib = "CP",
    filterThreshold = 0.5
)

# Group by time point
temporal_effects <- results %>%
    group_by(TargetTime) %>%
    summarize(
        n_drugs = n_distinct(Target),
        mean_similarity = mean(Similarity),
        .groups = "drop"
    ) %>%
    arrange(TargetTime)

cat("Time-dependent drug effects:\n")
print(temporal_effects)
```

### Case 4: Dose-Response Patterns

``` r

# Examine concentration-dependent effects
drug_doses <- results %>%
    filter(!is.na(TargetConcentration)) %>%
    group_by(Target, TargetConcentration) %>%
    summarize(
        mean_similarity = mean(Similarity),
        .groups = "drop"
    ) %>%
    group_by(Target) %>%
    filter(n() >= 2) %>% # Multiple concentrations tested
    ungroup()

# Plot dose-response for top drugs
top_drugs_dose <- drug_doses %>%
    group_by(Target) %>%
    slice_head(n = 5) %>%
    ungroup()

ggplot(top_drugs_dose, aes(x = TargetConcentration, y = mean_similarity)) +
    geom_line(aes(color = Target)) +
    geom_point(aes(color = Target), size = 2) +
    labs(
        title = "Dose-Response Patterns",
        x = "Concentration",
        y = "Mean Similarity"
    ) +
    theme_minimal()
```

## Integration with Experimental Data

### Validating Predictions

``` r

# After getting computational predictions
predictions <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5
)

# Select top candidates for validation
validation_candidates <- predictions %>%
    filter(Similarity < -0.4, pValue < 0.05) %>%
    arrange(Similarity) %>%
    head(10)

# Export for experimental testing
# write.csv(validation_candidates, "tp53_rescue_candidates.csv")
```

### Comparing with Literature

``` r

# Known TP53-interacting drugs from literature
literature_drugs <- c("nutlin-3", "PRIMA-1", "APR-246")

# Check if they appear in predictions
literature_hits <- predictions %>%
    filter(grepl(paste(literature_drugs, collapse = "|"),
        Target,
        ignore.case = TRUE
    ))

cat("Literature-validated drugs in predictions:\n")
print(literature_hits)
```

## Best Practices

### 1. Start Broad, Then Narrow

``` r

# Step 1: Broad discovery
broad_search <- investigateTarget(
    target = "EGFR",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5,
    similarityThreshold = 0.2
)

# Step 2: Focus on high-confidence
high_confidence <- broad_search %>%
    filter(abs(Similarity) > 0.4, pValue < 0.01)

# Step 3: Cell line-specific validation
validated <- investigateTarget(
    target = "EGFR",
    inputLib = "KD",
    outputLib = "CP",
    outputCellLines = c("A549", "MCF7"),
    filterThreshold = 0.6,
    similarityThreshold = 0.3
)
```

### 2. Consider Biological Context

- **Cell type**: Use appropriate cell line filters
- **Perturbation type**: KD vs OE have different implications
- **Magnitude**: Balance stringency with coverage
- **Direction**: Paired analysis for complex effects

### 3. Multiple Evidence Lines

``` r

# Evidence 1: KD → CP
evidence_kd <- investigateTarget(
    target = "BRCA1", inputLib = "KD", outputLib = "CP",
    filterThreshold = 0.5
)

# Evidence 2: OE → CP
evidence_oe <- investigateTarget(
    target = "BRCA1", inputLib = "OE", outputLib = "CP",
    filterThreshold = 0.5
)

# Find drugs appearing in both analyses
drugs_kd <- evidence_kd %>%
    filter(Similarity < -0.3) %>%
    pull(Target)
drugs_oe <- evidence_oe %>%
    filter(Similarity < -0.3) %>%
    pull(Target)

robust_candidates <- intersect(drugs_kd, drugs_oe)
cat("Drugs validated across KD and OE:", length(robust_candidates), "\n")
```

### 4. Document Parameters

Always record your analysis parameters:

``` r

analysis_params <- list(
    target = "TP53",
    input_library = "KD",
    output_library = "CP",
    filter_threshold = 0.5,
    similarity_threshold = 0.3,
    paired = TRUE,
    date = Sys.Date(),
    drugfindR_version = packageVersion("drugfindR")
)

# Save with results
# saveRDS(list(results = tp53_results, params = analysis_params),
#         "tp53_analysis.rds")
```

## Troubleshooting

### No Signatures Found

``` r

# If no signatures found for your target:

# 1. Check spelling and case
# Note: drugfindR handles case-insensitivity, but verify gene name

# 2. Try alternate names
# Some genes have aliases (e.g., ERBB2 vs HER2)

# 3. Check metadata
kd_metadata <- drugfindR:::.loadMetadata("KD")
available_genes <- unique(kd_metadata$Source)

# Search for your gene
"TP53" %in% available_genes
grep("TP53", available_genes, value = TRUE, ignore.case = TRUE)
```

### Too Few Results

``` r

# If results are sparse:

# 1. Lower thresholds
relaxed <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.3, # More permissive
    similarityThreshold = 0.1
)

# 2. Try different libraries
# KD → KD instead of KD → CP
alternate <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "KD",
    filterThreshold = 0.5
)

# 3. Remove cell line restrictions
unrestricted <- investigateTarget(
    target = "TP53",
    inputLib = "KD",
    outputLib = "CP",
    inputCellLines = NULL, # All cell lines
    outputCellLines = NULL,
    filterThreshold = 0.5
)
```

## Summary

Key takeaways for target investigation:

1.  **Choose right libraries**: KD/OE for genes, CP for drugs
2.  **Direction matters**: Positive = similar, negative = opposite
3.  **Use paired analysis**: For nuanced biological interpretation
4.  **Consider cell context**: Filter by relevant cell lines
5.  **Validate predictions**: Multiple evidence lines + experiments
6.  **Document thoroughly**: Record all parameters

## Next Steps

- **[Drug
  Repurposing](https://cogdisreslab.github.io/drugfindR/articles/drug-repurposing.md)**:
  Apply to disease signatures
- **[Getting
  Started](https://cogdisreslab.github.io/drugfindR/articles/getting-started.md)**:
  Quick reference guide
- **[Function
  Reference](https://cogdisreslab.github.io/drugfindR/reference/)**:
  Complete API

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
#> [1] purrr_1.2.2      tidyr_1.3.2      ggplot2_4.0.3    dplyr_1.2.1     
#> [5] drugfindR_1.0.0  BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] rappdirs_0.3.4      sass_0.4.10         generics_0.1.4     
#>  [4] stringi_1.8.7       hms_1.1.4           digest_0.6.39      
#>  [7] magrittr_2.0.5      RColorBrewer_1.1-3  evaluate_1.0.5     
#> [10] grid_4.6.1          bookdown_0.47       fastmap_1.2.0      
#> [13] jsonlite_2.0.0      DFplyr_1.6.0        BiocManager_1.30.27
#> [16] scales_1.4.0        httr2_1.2.3         textshaping_1.0.5  
#> [19] jquerylib_0.1.4     cli_3.6.6           rlang_1.3.0        
#> [22] withr_3.0.3         cachem_1.1.0        yaml_2.3.12        
#> [25] otel_0.2.0          tools_4.6.1         tzdb_0.5.0         
#> [28] BiocGenerics_0.58.1 curl_7.1.0          vctrs_0.7.3        
#> [31] R6_2.6.1            stats4_4.6.1        lifecycle_1.0.5    
#> [34] stringr_1.6.0       S4Vectors_0.50.1    fs_2.1.0           
#> [37] htmlwidgets_1.6.4   ragg_1.5.2          pkgconfig_2.0.3    
#> [40] desc_1.4.3          pkgdown_2.2.1       pillar_1.11.1      
#> [43] bslib_0.11.0        gtable_0.3.6        glue_1.8.1         
#> [46] systemfonts_1.3.2   xfun_0.60           tibble_3.3.1       
#> [49] tidyselect_1.2.1    knitr_1.51          farver_2.1.2       
#> [52] htmltools_0.5.9     rmarkdown_2.31      readr_2.2.0        
#> [55] compiler_4.6.1      S7_0.2.2
```

## References

1.  Subramanian A, et al. (2017). A Next Generation Connectivity Map:
    L1000 Platform and the First 1,000,000 Profiles. *Cell*,
    171(6):1437-1452.

2.  Lamb J, et al. (2006). The Connectivity Map: using gene-expression
    signatures to connect small molecules, genes, and disease.
    *Science*, 313(5795):1929-1935.

3.  Corsello SM, et al. (2020). Discovering the anticancer potential of
    non-oncology drugs by systematic viability profiling. *Nature
    Cancer*, 1(2):235-248.
