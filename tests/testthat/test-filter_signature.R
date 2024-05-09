# Test Input Validation

test_that("Not specifying threshold and prop causes error", {
    expect_error(filterSignature(exampleSignature()))
})

test_that("Specifying both threshold and prop causes error", {
    expect_error(filterSignature(exampleSignature(),
        threshold = 0.1, prop = 0.1
    ))
})

test_that("Impty signature when filtered is empty", {
    expect_identical(
        nrow(filterSignature(emptySignature(), threshold = 0.0)),
        0L
    )
})

test_that("Invalid signature direction causes error", {
    expect_error(filterSignature(exampleSignature(), direction = "invalid"))
})

test_that("More than two threshold values causes error", {
    expect_error(
        filterSignature(exampleSignature(), threshold = c(0.0, 0.1, 0.2))
    )
})

# Testing Filtering by Threshold

## Testing with one threshold value

test_that("Filter by threshold works", {
    expect_identical(
        nrow(filterSignature(exampleSignature(), threshold = 0.0)),
        978L
    )
})

test_that("Filter by Threshold works with non-zero threshold", {
    filtered <- filterSignature(exampleSignature(), threshold = 1.0)
    expect_true(all(abs(filtered[["Value_LogDiffExp"]]) >= 1.0))
})

test_that("Filter by Threshold works with non-zero threshold and explicity any", { # nolint: line_length_linter.
    filtered <- filterSignature(
        exampleSignature(),
        threshold = 1.0, direction = "any"
    )
    expect_true(all(
        abs(filtered[["Value_LogDiffExp"]]) >= 1.0
    ))
})

test_that("Filter by Threshold works on up-regulated genes", {
    filtered <- filterSignature(
        exampleSignature(),
        threshold = 1.0, direction = "up"
    )
    expect_true(all(filtered[["Value_LogDiffExp"]] >= 1.0))
})

test_that("Filter by Threshold works on down-regulated genes", {
    filtered <- filterSignature(
        exampleSignature(),
        threshold = 1.0, direction = "down"
    )
    expect_true(all(filtered[["Value_LogDiffExp"]] <= -1.0))
})

## Testing with two threshold values

test_that("Filter by Threshold works with two threshold values", {
    filtered <- filterSignature(exampleSignature(), threshold = c(-0.75, 1.0))
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filteredUp[["Value_LogDiffExp"]] >= 1.0))
    expect_true(all(filteredDown[["Value_LogDiffExp"]] <= -0.75))
    expect_identical(nrow(filtered), 458L)
})

test_that("Filter by Threshold works with two threshold values with explicit any", { # nolint: line_length_linter.
    filtered <- filterSignature(
        exampleSignature(),
        threshold = c(-0.75, 1.0), direction = "any"
    )
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filteredUp[["Value_LogDiffExp"]] >= 1.0))
    expect_true(all(filteredDown[["Value_LogDiffExp"]] <= -0.75))
    expect_identical(nrow(filtered), 458L)
})

test_that("Filter by Threshold works with two threshold values on up-regulated genes", { # nolint: line_length_linter.
    filtered <- filterSignature(
        exampleSignature(),
        threshold = c(-0.75, 1.0), direction = "up"
    )
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filteredUp[["Value_LogDiffExp"]] >= 1.0))
    expect_identical(nrow(filteredDown), 0L)
    expect_identical(nrow(filtered), 202L)
})

test_that("Filter by Threshold works with two threshold values on down-regulated genes", { # nolint: line_length_linter.
    filtered <- filterSignature(
        exampleSignature(),
        threshold = c(-0.75, 1.0), direction = "down"
    )
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_identical(nrow(filteredUp), 0L)
    expect_true(all(filteredDown[["Value_LogDiffExp"]] <= -0.75))
    expect_identical(nrow(filtered), 256L)
})

# Testing Filtering by Proportion

test_that("Filter by proportion works", {
    expect_identical(
        nrow(filterSignature(exampleSignature(), prop = 1.0)),
        978L
    )
})

test_that("Filter by proportion works with non-zero proportion", {
    exampleSignatueData <- exampleSignature()
    exampleUpThreshold <- quantile(
        exampleSignatueData[["Value_LogDiffExp"]], 0.9
    )
    exampleDownThreshold <- quantile(
        exampleSignatueData[["Value_LogDiffExp"]], 0.1
    )
    filtered <- filterSignature(exampleSignature(), prop = 0.1)
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filteredUp[["Value_LogDiffExp"]] >= exampleUpThreshold))
    expect_true(
        all(filteredDown[["Value_LogDiffExp"]] <= exampleDownThreshold)
    )
})

test_that("Filter by proportion works with non-zero proportion and explicit any", { # nolint: line_length_linter.
    exampleSignatueData <- exampleSignature()
    exampleUpThreshold <- quantile(
        exampleSignatueData[["Value_LogDiffExp"]], 0.9
    )
    exampleDownThreshold <- quantile(
        exampleSignatueData[["Value_LogDiffExp"]], 0.1
    )
    filtered <- filterSignature(
        exampleSignature(),
        prop = 0.1, direction = "any"
    )
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(
        all(filteredUp[["Value_LogDiffExp"]] >= exampleUpThreshold)
    )
    expect_true(
        all(filteredDown[["Value_LogDiffExp"]] <= exampleDownThreshold)
    )
})

test_that("Filter by proportion works on up-regulated genes", {
    exampleSignatueData <- exampleSignature()
    exampleUpThreshold <- quantile(
        exampleSignatueData[["Value_LogDiffExp"]], 0.9
    )
    filtered <- filterSignature(
        exampleSignature(),
        prop = 0.1, direction = "up"
    )
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filteredUp[["Value_LogDiffExp"]] >= exampleUpThreshold))
    expect_identical(nrow(filteredDown), 0L)
})

test_that("Filter by proportion works on down-regulated genes", {
    exampleSignatueData <- exampleSignature()
    exampleDownThreshold <- quantile(
        exampleSignatueData[["Value_LogDiffExp"]], 0.1
    )
    filtered <- filterSignature(
        exampleSignature(),
        prop = 0.1, direction = "down"
    )
    filteredUp <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filteredDown <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_identical(nrow(filteredUp), 0L)
    expect_true(
        all(filteredDown[["Value_LogDiffExp"]] <= exampleDownThreshold)
    )
})
