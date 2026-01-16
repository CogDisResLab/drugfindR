# Validate consensusConcordants input parameters

This internal function validates all input parameters for the
consensusConcordants function to ensure they meet the required
constraints.

## Usage

``` r
.validateConsensusConcordantsInput(dots, paired, cutoff, cellLine)
```

## Arguments

- dots:

  A list of dataframes passed via ... parameter.

- paired:

  Logical indicating whether paired analysis is requested.

- cutoff:

  Numeric similarity cutoff value.

- cellLine:

  Character vector of cell lines, or NULL.

## Value

Invisible NULL. The function throws an error if validation fails.

## Details

This function performs the following validations:

1.  Ensures paired analysis has exactly two dataframes

2.  Ensures unpaired analysis has exactly one dataframe

3.  Validates cutoff is numeric and within reasonable range

4.  Validates cellLine parameter format

## Examples

``` r
# Valid calls (no errors)
testData <- data.frame(similarity = c(0.5, -0.3), compound = c("A", "B"))
.validateConsensusConcordantsInput(list(testData), FALSE, 0.3, NULL)
#> Error in .validateConsensusConcordantsInput(list(testData), FALSE, 0.3,     NULL): could not find function ".validateConsensusConcordantsInput"
.validateConsensusConcordantsInput(list(testData, testData), TRUE, 0.3, "A375")
#> Error in .validateConsensusConcordantsInput(list(testData, testData),     TRUE, 0.3, "A375"): could not find function ".validateConsensusConcordantsInput"

# Invalid calls (will throw errors)
.validateConsensusConcordantsInput(list(), FALSE, 0.3, NULL) # No data
#> Error in .validateConsensusConcordantsInput(list(), FALSE, 0.3, NULL): could not find function ".validateConsensusConcordantsInput"
.validateConsensusConcordantsInput(list(testData), TRUE, 0.3, NULL) # Paired needs 2 dataframes
#> Error in .validateConsensusConcordantsInput(list(testData), TRUE, 0.3,     NULL): could not find function ".validateConsensusConcordantsInput"
```
