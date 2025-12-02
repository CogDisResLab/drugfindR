# Get concordant signatures from iLINCS database **\[stable\]**

This function queries the iLINCS (Integrative Library of Integrated
Network-based Cellular Signatures) database to find signatures that are
concordant (similar) to a given input signature.

## Usage

``` r
getConcordants(signature, ilincsLibrary = "CP")
```

## Arguments

- signature:

  A data.frame, tibble, or S4Vectors::DataFrame containing the signature
  data. Must conform to iLINCS signature structure with columns: \*
  `signatureID`: Signature identifier \* `ID_geneid`: Gene IDs \*
  `Name_GeneSymbol`: Gene symbols \* `Value_LogDiffExp`: Log fold-change
  values \* `Significance_pvalue`: Statistical significance p-values

  Use
  [`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.md)
  to ensure proper formatting.

- ilincsLibrary:

  Character string specifying the iLINCS library to search. Must be one
  of: \* `"CP"`: Chemical Perturbagen library (default) \* `"KD"`:
  Knockdown library \* `"OE"`: Overexpression library

## Value

A data structure containing concordant signatures. The return type
matches the input signature type: \* tibble for data.frame or tibble
inputs \* S4Vectors::DataFrame for DataFrame inputs

Contains the following columns: \* `signatureid`: Unique signature
identifier \* `compound` or `treatment`: Drug/treatment name \*
`concentration`: Drug concentration (CP library only) \* `time`:
Treatment duration \* `cellline`: Cell line used \* `similarity`:
Similarity score (rounded to 8 decimal places) \* `pValue`: Statistical
significance (rounded to 20 decimal places) \* `sig_direction`:
Signature direction ("Up", "Down", or "Any")

## Details

The function performs the following steps:

1.  Validates input parameters

2.  Creates a temporary file with signature data

3.  Detects signature direction from expression values

4.  Sends a multipart POST request to the iLINCS API

5.  Processes the JSON response into a standardized tibble

6.  Cleans up temporary files

The signature direction is determined as follows:

- `"Up"`: All expression values are greater than or equal to zero

- `"Down"`: All expression values are less than or equal to zero

- `"Any"`: Mixed positive and negative values

## API Details

This function interfaces with the iLINCS web service API. The signature
is uploaded as a tab-separated file and analyzed against the specified
library. Results are returned as JSON and parsed into a tibble.

## Error Handling

The function will stop execution with informative error messages for:

- Invalid signature data types (must be data.frame, tibble, or
  DataFrame)

- Invalid signature structure (missing required columns, wrong order,
  etc.)

- Missing values in signature data

- Unsupported iLINCS library names

- HTTP errors from the iLINCS API

- Invalid or empty API responses

## References

iLINCS Portal: <http://www.ilincs.org/>

Pilarczyk et al. (2020). Connecting omics signatures and revealing
biological mechanisms with iLINCS. Nature Communications, 11(1), 4058.

## See also

`[ prepareSignature() ]` for signature preparation,
`[ filterSignature() ]` for signature filtering,
`[ investigateSignature() ]` for signature investigation

## Examples

``` r
# Input validation examples (no API calls)
# These demonstrate proper signature structure
mockSig <- data.frame(
    signatureID = rep("TEST", 3),
    ID_geneid = c("123", "456", "789"),
    Name_GeneSymbol = c("TP53", "MYC", "EGFR"),
    Value_LogDiffExp = c(1.5, -2.0, 0.8)
)

# Validate library parameter (should produce error)
tryCatch(
    getConcordants(mockSig, ilincsLibrary = "INVALID"),
    error = function(e) message("Expected error: invalid library")
)
#> Expected error: invalid library

# This example requires network access to the iLINCS API

# Load example differential expression data
dge_file <- system.file("extdata", "dCovid_diffexp.tsv",
    package = "drugfindR"
)
dge_data <- read.delim(dge_file)

# Prepare signature to ensure proper structure
signature <- prepareSignature(
    dge_data[1:50, ],
    geneColumn = "hgnc_symbol",
    logfcColumn = "logFC",
    pvalColumn = "PValue"
)

# Find concordant chemical perturbagens
cpConcordants <- getConcordants(signature, ilincsLibrary = "CP")
head(cpConcordants)
#> # A tibble: 6 × 9
#>   signatureid    treatment       concentration time  cellline similarity  pValue
#>   <chr>          <chr>           <chr>         <chr> <chr>         <dbl>   <dbl>
#> 1 LINCSCP_83965  SCHEMBL15556278 0.37uM        24h   ASC.C         1.000 2.18e-6
#> 2 LINCSCP_59213  CHEMBL2361430   10uM          24h   VCAP         -1.000 3.91e-6
#> 3 LINCSCP_6945   Coleoforsin     10uM          6h    A549         -1.000 8.70e-6
#> 4 LINCSCP_251399 Tozasertib      1uM           24h   MCF7         -1.000 1.48e-5
#> 5 LINCSCP_26754  XMD-16144       10uM          6h    HT29         -1.000 1.72e-5
#> 6 LINCSCP_46851  AC1M43RN        10uM          24h   PC3           1.000 1.91e-5
#> # ℹ 2 more variables: sig_direction <chr>, sig_type <chr>

# Find concordant knockdown signatures
kdConcordants <- getConcordants(signature, ilincsLibrary = "KD")
head(kdConcordants)
#> # A tibble: 6 × 9
#>   signatureid   treatment concentration time  cellline similarity    pValue
#>   <chr>         <chr>     <chr>         <chr> <chr>         <dbl>     <dbl>
#> 1 LINCSKD_20320 RNF145    NA            96 h  HT29         -1.000 0.0000193
#> 2 LINCSKD_17785 ZNF267    NA            96 h  HEPG2        -1.000 0.0000566
#> 3 LINCSKD_12636 HNF4G     NA            96 h  HCC515        1.000 0.0000624
#> 4 LINCSKD_14114 ZNF621    NA            96 h  HCC515        1.000 0.0000788
#> 5 LINCSKD_6737  TAOK3     NA            96 h  A549         -1.000 0.000107 
#> 6 LINCSKD_2259  PCSK7     NA            96 h  A375         -1.000 0.000110 
#> # ℹ 2 more variables: sig_direction <chr>, sig_type <chr>

# Find concordant overexpression signatures
oeConcordants <- getConcordants(signature, ilincsLibrary = "OE")
head(oeConcordants)
#> # A tibble: 6 × 9
#>   signatureid   treatment concentration time  cellline similarity    pValue
#>   <chr>         <chr>     <chr>         <chr> <chr>         <dbl>     <dbl>
#> 1 LINCSOE_6126  GNB1      NA            96 h  HA1E          1.000 0.0000234
#> 2 LINCSOE_4496  INHBE     NA            96 h  A549         -1.000 0.000329 
#> 3 LINCSOE_5490  GPR65     NA            96 h  HA1E         -1.000 0.000450 
#> 4 LINCSOE_10043 SMAD3     NA            48 h  HEK293T      -1.000 0.000500 
#> 5 LINCSOE_10325 ZRSR1     NA            48 h  HEK293T      -0.999 0.000510 
#> 6 LINCSOE_216   HACD1     NA            96 h  A375          0.999 0.000517 
#> # ℹ 2 more variables: sig_direction <chr>, sig_type <chr>

# Works with different data frame types
signatureDf <- as.data.frame(signature)
cpConcordantsDf <- getConcordants(signatureDf, "CP")

# Works with S4Vectors::DataFrame
signatureDataFrame <- S4Vectors::DataFrame(signature)
cpConcordantsDataFrame <- getConcordants(signatureDataFrame, "CP")
# Returns S4Vectors::DataFrame to match input type
```
