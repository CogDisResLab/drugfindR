# Group concordants by target and select maximum similarity entries

This internal function groups concordants data by target (compound or
treatment) and retains only the entries with maximum absolute similarity
for each target.

## Usage

``` r
.groupByTargetAndSelectMax(concordants)
```

## Arguments

- concordants:

  A dataframe containing filtered concordants data.

## Value

A dataframe with deduplicated targets, keeping maximum similarity
entries.

## Details

This function:

1.  Groups by treatment or compound columns (whichever is available)

2.  For each group, retains only entries with maximum absolute
    similarity

3.  Handles ties by keeping all tied entries

4.  Preserves the structure for downstream processing

## Examples

``` r
testData <- data.frame(
    compound = c("A", "A", "B", "B"),
    similarity = c(0.5, 0.8, -0.3, -0.7),
    cellline = c("A375", "PC3", "A375", "PC3")
)
grouped <- .groupByTargetAndSelectMax(testData)
#> Error in .groupByTargetAndSelectMax(testData): could not find function ".groupByTargetAndSelectMax"
# Returns entries with max |similarity| for each compound
```
