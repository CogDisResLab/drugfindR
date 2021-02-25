library(readr)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(stringr)

concordant_col_spec <- cols(.default = col_character(),
                            similarity = col_double(),
                            pValue = col_double())

consensus_col_spec <- cols(
  .default = col_character(),
  Similarity = col_double()
)

unpaired <- list.files("reference/", "concordant.*any") %>%
  file.path("reference/", .)

unpaired_data <-
  map(unpaired, ~ read_csv(.x, col_types = concordant_col_spec))

up_paired <- list.files("reference/", "concordant.*up") %>%
  file.path("reference/", .)

up_paired_data <-
  map(up_paired, ~ read_csv(.x, col_types = concordant_col_spec))

down_paired <- list.files("reference/", "concordant.*down") %>%
  file.path("reference/", .)

down_paired_data <-
  map(down_paired, ~ read_csv(.x, col_types = concordant_col_spec))

paired_data <- map2(up_paired_data, down_paired_data, list)

unpaired_pm_files <-
  expand_grid(
    unpaired,
    paired = FALSE,
    cutoff = c(0.2, 0.4),
    cell_line = list(NULL, "HA1E", c("MCF7, HA1E")),
    discordant = c(TRUE, FALSE)
  ) %>% mutate(
    paired = "unpaired",
    cell_line = case_when(
      cell_line == "NULL" ~ "null",
      cell_line == "HA1E" ~ "single",
      TRUE ~ "multiple"
    ),
    consensus = "consensus",
    discordant = if_else(discordant, "discordant", "concordant"),
    source = str_extract(unpaired, "(cp|oe|kd)_signature"),
    target = str_extract(unpaired, "(CP|OE|KD)")
  ) %>% select(consensus, source, target, paired, cutoff, cell_line, discordant) %>%
  unite(filename) %>%
  pull(filename) %>%
  str_c(., "csv", sep = ".") %>%
  file.path("reference", .)

paired_pm_files <-
  expand_grid(
    up_paired,
    paired = TRUE,
    cutoff = c(0.2, 0.4),
    cell_line = list(NULL, "HA1E", c("MCF7, HA1E")),
    discordant = c(TRUE, FALSE)
  ) %>% mutate(
    paired = "paired",
    cell_line = case_when(
      cell_line == "NULL" ~ "null",
      cell_line == "HA1E" ~ "single",
      TRUE ~ "multiple"
    ),
    consensus = "consensus",
    discordant = if_else(discordant, "discordant", "concordant"),
    source = str_extract(up_paired, "(cp|oe|kd)_signature"),
    target = str_extract(up_paired, "(CP|OE|KD)")
  ) %>%
  select(consensus, source, target, paired, cutoff, cell_line, discordant) %>%
  unite(filename) %>%
  pull(filename) %>%
  str_c(., "csv", sep = ".") %>%
  file.path("reference", .)

# Testing unpaired data

# Testing CP-CP Comparison for 0.2

test_that("CP-CP consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[1], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[2], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[3], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[4], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[5], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[6], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-CP Comparison for 0.4

test_that("CP-CP consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[7], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[8], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[9], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[10], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[11], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[12], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[1], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-KD Comparison for 0.2

test_that("CP-KD consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[13], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[14], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[15], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[16], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[17], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[18], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-KD Comparison for 0.4

test_that("CP-KD consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[19], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[20], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[21], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[22], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[23], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[24], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[2], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-OE Comparison for 0.2

test_that("CP-OE consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[25], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[26], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[27], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[28], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[29], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[30], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-OE Comparison for 0.4

test_that("CP-OE consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[31], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[32], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[33], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[34], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[35], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[36], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[3], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-CP Comparison for 0.2

test_that("KD-CP consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[37], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[38], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[39], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[40], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[41], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[42], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-CP Comparison for 0.4

test_that("KD-CP consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[43], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[44], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[45], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[46], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[47], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[48], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[4], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-KD Comparison for 0.2

test_that("KD-KD consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[49], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[50], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[51], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[52], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[53], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[54], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-KD Comparison for 0.4

test_that("KD-KD consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[55], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[56], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[57], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[58], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[59], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[60], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[5], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-OE Comparison for 0.2

test_that("KD-OE consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[61], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[62], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[63], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[64], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[65], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[66], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-OE Comparison for 0.4

test_that("KD-OE consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[67], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[68], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[69], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[70], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[71], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[72], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[6], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-CP Comparison for 0.2

test_that("OE-CP consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[73], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[74], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[75], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[76], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[77], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[78], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-CP Comparison for 0.4

test_that("OE-CP consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[79], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[80], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[81], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[82], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[83], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[84], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[7], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-KD Comparison for 0.2

test_that("OE-KD consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[85], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[86], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[87], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[88], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[89], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[90], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-KD Comparison for 0.4

test_that("OE-KD consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[91], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[92], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[93], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[94], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[95], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[96], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[8], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-OE Comparison for 0.2

test_that("OE-OE consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[97], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[98], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[99], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[100], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[101], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(unpaired_pm_files[102], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-OE Comparison for 0.4

test_that("OE-OE consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[103], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[104], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[105], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[106], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[107], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(unpaired_pm_files[108], col_types = consensus_col_spec)
  example <- consensus_concordants(unpaired_data[9], paired = FALSE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})


# Testing paired data

# Testing CP-CP Comparison for 0.2

test_that("CP-CP consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[1], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[2], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[3], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[4], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[5], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[6], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-CP Comparison for 0.4

test_that("CP-CP consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[7], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[8], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[9], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[10], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[11], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-CP consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[12], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[1]][[1]], paired_data[[1]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-KD Comparison for 0.2

test_that("CP-KD consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[13], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[14], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[15], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[16], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[17], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[18], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-KD Comparison for 0.4

test_that("CP-KD consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[19], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[20], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[21], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[22], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[23], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-KD consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[24], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[2]][[1]], paired_data[[2]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-OE Comparison for 0.2

test_that("CP-OE consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[25], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[26], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[27], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[28], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[29], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[30], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing CP-OE Comparison for 0.4

test_that("CP-OE consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[31], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[32], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[33], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[34], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[35], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("CP-OE consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[36], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[3]][[1]], paired_data[[3]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-CP Comparison for 0.2

test_that("KD-CP consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[37], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[38], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[39], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[40], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[41], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[42], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-CP Comparison for 0.4

test_that("KD-CP consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[43], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[44], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[45], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[46], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[47], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-CP consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[48], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[4]][[1]], paired_data[[4]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-KD Comparison for 0.2

test_that("KD-KD consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[49], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[50], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[51], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[52], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[53], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[54], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-KD Comparison for 0.4

test_that("KD-KD consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[55], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[56], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[57], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[58], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[59], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-KD consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[60], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[5]][[1]], paired_data[[5]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-OE Comparison for 0.2

test_that("KD-OE consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[61], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[62], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[63], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[64], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[65], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[66], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing KD-OE Comparison for 0.4

test_that("KD-OE consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[67], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[68], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[69], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[70], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[71], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("KD-OE consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[72], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[6]][[1]], paired_data[[6]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-CP Comparison for 0.2

test_that("OE-CP consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[73], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[74], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[75], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[76], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[77], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[78], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-CP Comparison for 0.4

test_that("OE-CP consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[79], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[80], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[81], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[82], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[83], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-CP consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[84], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[7]][[1]], paired_data[[7]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-KD Comparison for 0.2

test_that("OE-KD consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[85], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[86], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[87], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[88], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[89], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[90], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-KD Comparison for 0.4

test_that("OE-KD consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[91], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[92], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[93], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[94], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[95], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-KD consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[96], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[8]][[1]], paired_data[[8]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-OE Comparison for 0.2

test_that("OE-OE consensus discordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[97], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for no cell line and 0.2", {
  reference <- read_csv(paired_pm_files[98], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.2, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[99], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for one cell line and 0.2", {
  reference <- read_csv(paired_pm_files[100], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[101], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.2))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for multiple cell line and 0.2", {
  reference <- read_csv(paired_pm_files[102], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.2, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.2))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

# Testing OE-OE Comparison for 0.4

test_that("OE-OE consensus discordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[103], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = TRUE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for no cell line and 0.4", {
  reference <- read_csv(paired_pm_files[104], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.4, cell_line = NULL, discordant = FALSE)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[105], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for one cell line and 0.4", {
  reference <- read_csv(paired_pm_files[106], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 1 && cell_line == "HA1E") || length(cell_line) == 0)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus discordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[107], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = TRUE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity < -0.4))
  expect_true(all(example$Similarity < 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

test_that("OE-OE consensus concordants are okay for multiple cell line and 0.4", {
  reference <- read_csv(paired_pm_files[108], col_types = consensus_col_spec)
  example <- consensus_concordants(paired_data[[9]][[1]], paired_data[[9]][[2]], paired = TRUE, cutoff = 0.4, cell_line = c("HA1E", "MCF7"), discordant = FALSE)
  cell_line <- unique(example$TargetCellLine)
  expect_equal(nrow(example), nrow(reference))
  expect_true(all(example$Similarity > 0.4))
  expect_true(all(example$Similarity > 0))
  expect_true((length(cell_line) == 2 && cell_line %in% c("HA1E", "MCF7")) || length(cell_line) < 2)
  expect_equal(example, reference, ignore_attr = TRUE)
})

