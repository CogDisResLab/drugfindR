# Test Input Validation

test_that("Not specifying threshold and prop causes error", {
    expect_error(filter_signature(example_signature()))
})

test_that("Specifying both threshold and prop causes error", {
    expect_error(filter_signature(kd_signature, threshold = 0.1, prop = 0.1))
})

test_that("Impty signature when filtered is empty", {
    expect_identical(
        nrow(filter_signature(empty_signature(), threshold = 0.0)),
        0L
    )
})

test_that("Invalid signature direction causes error", {
    expect_error(filter_signature(example_signature(), direction = "invalid"))
})

test_that("More than two threshold values causes error", {
    expect_error(
        filter_signature(example_signature(), threshold = c(0.0, 0.1, 0.2))
    )
})

# Testing Filtering by Threshold

## Testing with one threshold value

test_that("Filter by threshold works", {
    expect_identical(
        nrow(filter_signature(example_signature(), threshold = 0.0)),
        978L
    )
})

test_that("Filter by Threshold works with non-zero threshold", {
    filtered <- filter_signature(example_signature(), threshold = 1.0)
    expect_true(all(abs(filtered[["Value_LogDiffExp"]]) >= 1.0))
})

test_that("Filter by Threshold works with non-zero threshold and explicity any", { # nolint: line_length_linter.
    filtered <- filter_signature(
        example_signature(),
        threshold = 1.0, direction = "any"
    )
    expect_true(all(
        abs(filtered[["Value_LogDiffExp"]]) >= 1.0
    ))
})

test_that("Filter by Threshold works on up-regulated genes", {
    filtered <- filter_signature(
        example_signature(),
        threshold = 1.0, direction = "up"
    )
    expect_true(all(filtered[["Value_LogDiffExp"]] >= 1.0))
})

test_that("Filter by Threshold works on down-regulated genes", {
    filtered <- filter_signature(
        example_signature(),
        threshold = 1.0, direction = "down"
    )
    expect_true(all(filtered[["Value_LogDiffExp"]] <= -1.0))
})

## Testing with two threshold values

test_that("Filter by Threshold works with two threshold values", {
    filtered <- filter_signature(example_signature(), threshold = c(-0.75, 1.0))
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filtered_up[["Value_LogDiffExp"]] >= 1.0))
    expect_true(all(filtered_down[["Value_LogDiffExp"]] <= -0.75))
    expect_identical(nrow(filtered), 458L)
})

test_that("Filter by Threshold works with two threshold values with explicit any", { # nolint: line_length_linter.
    filtered <- filter_signature(
        example_signature(),
        threshold = c(-0.75, 1.0), direction = "any"
    )
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filtered_up[["Value_LogDiffExp"]] >= 1.0))
    expect_true(all(filtered_down[["Value_LogDiffExp"]] <= -0.75))
    expect_identical(nrow(filtered), 458L)
})

test_that("Filter by Threshold works with two threshold values on up-regulated genes", { # nolint: line_length_linter.
    filtered <- filter_signature(
        example_signature(),
        threshold = c(-0.75, 1.0), direction = "up"
    )
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filtered_up[["Value_LogDiffExp"]] >= 1.0))
    expect_identical(nrow(filtered_down), 0L)
    expect_identical(nrow(filtered), 202L)
})

test_that("Filter by Threshold works with two threshold values on down-regulated genes", { # nolint: line_length_linter.
    filtered <- filter_signature(
        example_signature(),
        threshold = c(-0.75, 1.0), direction = "down"
    )
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_identical(nrow(filtered_up), 0L)
    expect_true(all(filtered_down[["Value_LogDiffExp"]] <= -0.75))
    expect_identical(nrow(filtered), 256L)
})

# Testing Filtering by Proportion

test_that("Filter by proportion works", {
    expect_identical(
        nrow(filter_signature(example_signature(), prop = 1.0)),
        978L
    )
})

test_that("Filter by proportion works with non-zero proportion", {
    example_signatue_data <- example_signature()
    example_up_threshold <- quantile(
        example_signatue_data[["Value_LogDiffExp"]], 0.9
    )
    example_down_threshold <- quantile(
        example_signatue_data[["Value_LogDiffExp"]], 0.1
    )
    filtered <- filter_signature(example_signature(), prop = 0.1)
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filtered_up[["Value_LogDiffExp"]] >= example_up_threshold))
    expect_true(
        all(filtered_down[["Value_LogDiffExp"]] <= example_down_threshold)
    )
})

test_that("Filter by proportion works with non-zero proportion and explicit any", { # nolint: line_length_linter.
    example_signatue_data <- example_signature()
    example_up_threshold <- quantile(
        example_signatue_data[["Value_LogDiffExp"]], 0.9
    )
    example_down_threshold <- quantile(
        example_signatue_data[["Value_LogDiffExp"]], 0.1
    )
    filtered <- filter_signature(
        example_signature(),
        prop = 0.1, direction = "any"
    )
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(
        all(filtered_up[["Value_LogDiffExp"]] >= example_up_threshold)
    )
    expect_true(
        all(filtered_down[["Value_LogDiffExp"]] <= example_down_threshold)
    )
})

test_that("Filter by proportion works on up-regulated genes", {
    example_signatue_data <- example_signature()
    example_up_threshold <- quantile(
        example_signatue_data[["Value_LogDiffExp"]], 0.9
    )
    filtered <- filter_signature(
        example_signature(),
        prop = 0.1, direction = "up"
    )
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_true(all(filtered_up[["Value_LogDiffExp"]] >= example_up_threshold))
    expect_identical(nrow(filtered_down), 0L)
})

test_that("Filter by proportion works on down-regulated genes", {
    example_signatue_data <- example_signature()
    example_down_threshold <- quantile(
        example_signatue_data[["Value_LogDiffExp"]], 0.1
    )
    filtered <- filter_signature(
        example_signature(),
        prop = 0.1, direction = "down"
    )
    filtered_up <- filtered[filtered[["Value_LogDiffExp"]] >= 0.0, ]
    filtered_down <- filtered[filtered[["Value_LogDiffExp"]] <= 0.0, ]
    expect_identical(nrow(filtered_up), 0L)
    expect_true(
        all(filtered_down[["Value_LogDiffExp"]] <= example_down_threshold)
    )
})
