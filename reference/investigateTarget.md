# Investigate concordant signatures for a gene or drug **\[stable\]**

Given the name of a target (gene knockdown/overexpression or compound)
this high–level convenience wrapper:

## Usage

``` r
investigateTarget(
  target,
  inputLib,
  outputLib,
  filterThreshold = 0.85,
  similarityThreshold = 0.321,
  paired = TRUE,
  inputCellLines = NULL,
  outputCellLines = NULL
)
```

## Arguments

- target:

  Character scalar. Gene symbol (for KD/OE libraries), or drug /
  compound name (for CP library) used to locate source signatures.

- inputLib:

  Character ("OE", "KD", or "CP"). Library from which source signatures
  for `target` are drawn.

- outputLib:

  Character ("OE", "KD", or "CP"). Library queried for concordant
  signatures.

- filterThreshold:

  Numeric in \\0,1\]. Minimum absolute (or directional) change used to
  retain genes in each source signature prior to concordance. Default
  `0.85` is conservative; consider lowering (e.g. `0.5`) for broader
  coverage.

- similarityThreshold:

  Numeric in \[0,1\]. Minimum similarity score retained in the final
  consensus result set. Default `0.321` (~ upper third).

- paired:

  Logical. If `TRUE` (default) computes concordance separately for up
  and down regulated gene sets; if `FALSE` uses all selected genes
  together.

- inputCellLines:

  Optional character vector restricting the search for source signatures
  to specified cell line(s). If `NULL` all available are considered.

- outputCellLines:

  Optional character vector restricting target signatures (during
  consensus formation) to specified cell line(s). If `NULL` all are
  considered.

## Value

A tibble (data frame) with one row per consensus concordant target
signature. Typical columns include:

- `Source` / `Target` – gene or compound names.

- `Similarity` – numeric concordance score in \[-1,1\].

- `SourceSignature`, `TargetSignature` – iLINCS signature identifiers.

- `SourceCellLine`, `TargetCellLine` – originating cell lines (if
  applicable).

- `SourceConcentration`, `TargetConcentration` – dosing information for
  CP.

- `SourceTime`, `TargetTime` – time point metadata.

## Details

1.  Locates iLINCS source signatures for the target in the specified
    input library.

2.  Optionally filters by source cell line(s).

3.  Retrieves each source signature and filters genes by direction and
    magnitude.

4.  Queries iLINCS for concordant signatures in the chosen output
    library.

5.  Computes paired or unpaired consensus concordance across up/down
    regulated sets.

6.  Returns an augmented tibble of similarity scores and rich
    source/target metadata.

The paired workflow evaluates concordance separately for up‑ and
down‑regulated genes and then combines (via
[`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md))
the two result sets. When `paired = FALSE` a single aggregate signature
(direction = "any") is used.

Network access: This function performs remote API calls (unless tests
are run under a mocking context such as
[`httptest2::with_mock_api()`](https://enpiar.com/httptest2/reference/with_mock_api.html)).
Examples are wrapped in `\donttest{}` to avoid false negatives on CRAN /
Bioconductor builders without network access.

Errors are raised if:

- No source signatures match `target` in the requested `inputLib` (empty
  set).

- Invalid library codes are supplied.

Internally this function orchestrates:
[`getSignature()`](https://cogdisreslab.github.io/drugfindR/reference/getSignature.md),
[`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md),
[`getConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/getConcordants.md)
and
[`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md).
It returns a vertically concatenated result across all matching source
signatures.

## Thresholds

- `filterThreshold` controls gene selection within each source
  signature. It is passed to
  [`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md)
  as the absolute (or directional) threshold.

- `similarityThreshold` is applied when forming the consensus
  concordants to discard low similarity entries.

## See also

[`getSignature()`](https://cogdisreslab.github.io/drugfindR/reference/getSignature.md),
[`filterSignature()`](https://cogdisreslab.github.io/drugfindR/reference/filterSignature.md),
[`getConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/getConcordants.md),
[`consensusConcordants()`](https://cogdisreslab.github.io/drugfindR/reference/consensusConcordants.md),
[`prepareSignature()`](https://cogdisreslab.github.io/drugfindR/reference/prepareSignature.md)
for lower‑level operations.

## Examples

``` r
# Input validation examples (no API calls)
# Demonstrate library parameter validation
tryCatch(
    investigateTarget(target = "TP53", inputLib = "INVALID", outputLib = "CP"),
    error = function(e) message("Expected error: invalid inputLib")
)
#> Expected error: invalid inputLib

tryCatch(
    investigateTarget(target = "TP53", inputLib = "KD", outputLib = "INVALID"),
    error = function(e) message("Expected error: invalid outputLib")
)
#> Expected error: invalid outputLib

# \donttest{
# This function makes multiple API calls to iLINCS and may take several minutes
# Basic paired investigation of a knockdown signature against compound library
set.seed(1)
res <- investigateTarget(
    target = "AATK",
    inputLib = "KD",
    outputLib = "CP",
    filterThreshold = 0.5,
    similarityThreshold = 0.3,
    paired = TRUE
)
head(res)
#> # A tibble: 6 × 13
#>   Source Target        Similarity SourceSignature SourceCellLine SourceTime
#>   <chr>  <chr>              <dbl> <chr>           <chr>          <chr>     
#> 1 AATK   Withaferin A       0.572 LINCSKD_29946   VCAP           120 h     
#> 2 AATK   Adapalene          0.552 LINCSKD_29946   VCAP           120 h     
#> 3 AATK   BSPBIO_001592     -0.552 LINCSKD_29946   VCAP           120 h     
#> 4 AATK   BRD-K04759195      0.549 LINCSKD_29946   VCAP           120 h     
#> 5 AATK   Cabazitaxel       -0.545 LINCSKD_29946   VCAP           120 h     
#> 6 AATK   Erlotinib          0.543 LINCSKD_29946   VCAP           120 h     
#> # ℹ 7 more variables: TargetSignature <chr>, TargetCellLine <chr>,
#> #   TargetConcentration <chr>, TargetTime <chr>, InputSigDirection <chr>,
#> #   SignatureType <chr>, pValue <dbl>

# Unpaired (aggregate) workflow — often faster, returns a single consensus set
res_unpaired <- investigateTarget(
    target = "AATK", inputLib = "KD", outputLib = "CP",
    filterThreshold = 0.5, similarityThreshold = 0.3, paired = FALSE
)
head(res_unpaired)
#> # A tibble: 6 × 13
#>   Source Target          Similarity SourceSignature SourceCellLine SourceTime
#>   <chr>  <chr>                <dbl> <chr>           <chr>          <chr>     
#> 1 AATK   Spironolactone       0.671 LINCSKD_29946   VCAP           120 h     
#> 2 AATK   SCHEMBL16226502      0.635 LINCSKD_29946   VCAP           120 h     
#> 3 AATK   Danazol              0.632 LINCSKD_29946   VCAP           120 h     
#> 4 AATK   Simvastatin         -0.622 LINCSKD_29946   VCAP           120 h     
#> 5 AATK   ST056792             0.608 LINCSKD_29946   VCAP           120 h     
#> 6 AATK   BRD-K73434853        0.607 LINCSKD_29946   VCAP           120 h     
#> # ℹ 7 more variables: TargetSignature <chr>, TargetCellLine <chr>,
#> #   TargetConcentration <chr>, TargetTime <chr>, InputSigDirection <chr>,
#> #   SignatureType <chr>, pValue <dbl>

# Restrict source signatures to specific cell lines (if available)
# and target signatures to a subset of cell lines during consensus
res_filtered <- investigateTarget(
    target = "AATK", inputLib = "KD", outputLib = "CP",
    outputCellLines = c("MCF7"),
    filterThreshold = 0.5, similarityThreshold = 0.3
)
head(res_filtered)
#> # A tibble: 6 × 13
#>   Source Target       Similarity SourceSignature SourceCellLine SourceTime
#>   <chr>  <chr>             <dbl> <chr>           <chr>          <chr>     
#> 1 AATK   Erlotinib         0.543 LINCSKD_29946   VCAP           120 h     
#> 2 AATK   Wortmannin        0.503 LINCSKD_29946   VCAP           120 h     
#> 3 AATK   Dacarbazine       0.491 LINCSKD_29946   VCAP           120 h     
#> 4 AATK   Tadalafil         0.471 LINCSKD_29946   VCAP           120 h     
#> 5 AATK   MG-132            0.468 LINCSKD_29946   VCAP           120 h     
#> 6 AATK   Tacedinaline      0.468 LINCSKD_29946   VCAP           120 h     
#> # ℹ 7 more variables: TargetSignature <chr>, TargetCellLine <chr>,
#> #   TargetConcentration <chr>, TargetTime <chr>, InputSigDirection <chr>,
#> #   SignatureType <chr>, pValue <dbl>

# Using httptest2 (if installed) to mock network calls:
# httptest2::with_mock_api({
#   mock_res <- investigateTarget("AATK", "KD", "CP", filterThreshold = 0.5)
#   print(head(mock_res))
# })
# }
```
