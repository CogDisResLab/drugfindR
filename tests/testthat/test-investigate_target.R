library(readr)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(stringr)

s1 <- expand_grid(target = "Z36", input_lib = "CP", output_lib = c("KD", "OE", "CP"), discordant = c(FALSE, TRUE), similarity_threshold = c(0.2, 0.4), output_cell_lines = list(NULL, c("VCAP"), c("VCAP", "A375")))

s2 <- expand_grid(target = "ADRBK2", input_lib = "KD", output_lib = c("KD", "OE", "CP"), discordant = c(FALSE, TRUE), similarity_threshold = c(0.2, 0.4), output_cell_lines = list(NULL, c("VCAP"), c("VCAP", "A375")))

s3 <- expand_grid(target = "FBXO7", input_lib = "OE", output_lib = c("KD", "OE", "CP"), discordant = c(FALSE, TRUE), similarity_threshold = c(0.2, 0.4), output_cell_lines = list(NULL, c("VCAP"), c("VCAP", "A375")))

s <- bind_rows(s1, s2, s3)

ref_files <- s %>% mutate(discordant = if_else(discordant, "discordant", "concordant"), cell_lines = ifelse(output_cell_lines == "NULL", "null", NA), cell_lines = ifelse(output_cell_lines == "VCAP", "single", cell_lines), cell_lines = ifelse(is.na(cell_lines), "multiple", cell_lines)) %>% select(-output_cell_lines) %>% unite(filename) %>% pull(filename) %>% str_c("csv", sep = ".") %>% file.path("reference", .)

ref_data <- map(ref_files, read_csv)

# Testing CP-KD Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[1,], investigate_target)
  ref <- ref_data[[1]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[2,], investigate_target)
  ref <- ref_data[[2]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[3,], investigate_target)
  ref <- ref_data[[3]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[4,], investigate_target)
  ref <- ref_data[[4]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[5,], investigate_target)
  ref <- ref_data[[5]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[6,], investigate_target)
  ref <- ref_data[[6]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[7,], investigate_target)
  ref <- ref_data[[7]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[8,], investigate_target)
  ref <- ref_data[[8]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[9,], investigate_target)
  ref <- ref_data[[9]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[10,], investigate_target)
  ref <- ref_data[[10]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[11,], investigate_target)
  ref <- ref_data[[110]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[12,], investigate_target)
  ref <- ref_data[[12]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing CP-OE Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[13,], investigate_target)
  ref <- ref_data[[13]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[14,], investigate_target)
  ref <- ref_data[[14]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[15,], investigate_target)
  ref <- ref_data[[15]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[16,], investigate_target)
  ref <- ref_data[[16]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[17,], investigate_target)
  ref <- ref_data[[17]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[18,], investigate_target)
  ref <- ref_data[[18]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[19,], investigate_target)
  ref <- ref_data[[19]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[20,], investigate_target)
  ref <- ref_data[[20]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[21,], investigate_target)
  ref <- ref_data[[21]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[22,], investigate_target)
  ref <- ref_data[[22]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[23,], investigate_target)
  ref <- ref_data[[23]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[24,], investigate_target)
  ref <- ref_data[[24]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing CP-CP Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[25,], investigate_target)
  ref <- ref_data[[25]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[26,], investigate_target)
  ref <- ref_data[[26]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[27,], investigate_target)
  ref <- ref_data[[27]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[28,], investigate_target)
  ref <- ref_data[[28]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[29,], investigate_target)
  ref <- ref_data[[29]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[30,], investigate_target)
  ref <- ref_data[[30]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[31,], investigate_target)
  ref <- ref_data[[31]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[32,], investigate_target)
  ref <- ref_data[[32]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[33,], investigate_target)
  ref <- ref_data[[33]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[34,], investigate_target)
  ref <- ref_data[[34]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[35,], investigate_target)
  ref <- ref_data[[35]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[36,], investigate_target)
  ref <- ref_data[[36]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing KD-KD Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[37,], investigate_target)
  ref <- ref_data[[37]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[38,], investigate_target)
  ref <- ref_data[[38]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[39,], investigate_target)
  ref <- ref_data[[39]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[40,], investigate_target)
  ref <- ref_data[[40]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[41,], investigate_target)
  ref <- ref_data[[41]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[42,], investigate_target)
  ref <- ref_data[[42]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[43,], investigate_target)
  ref <- ref_data[[43]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[44,], investigate_target)
  ref <- ref_data[[44]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[45,], investigate_target)
  ref <- ref_data[[45]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[46,], investigate_target)
  ref <- ref_data[[46]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[47,], investigate_target)
  ref <- ref_data[[47]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[48,], investigate_target)
  ref <- ref_data[[48]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing KD-OE Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[49,], investigate_target)
  ref <- ref_data[[49]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[50,], investigate_target)
  ref <- ref_data[[50]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[51,], investigate_target)
  ref <- ref_data[[51]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[52,], investigate_target)
  ref <- ref_data[[52]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[53,], investigate_target)
  ref <- ref_data[[53]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[54,], investigate_target)
  ref <- ref_data[[54]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[55,], investigate_target)
  ref <- ref_data[[55]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[56,], investigate_target)
  ref <- ref_data[[56]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[57,], investigate_target)
  ref <- ref_data[[57]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[58,], investigate_target)
  ref <- ref_data[[58]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[59,], investigate_target)
  ref <- ref_data[[59]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[60,], investigate_target)
  ref <- ref_data[[60]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing KD-CP Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[61,], investigate_target)
  ref <- ref_data[[61]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[62,], investigate_target)
  ref <- ref_data[[62]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[63,], investigate_target)
  ref <- ref_data[[63]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[64,], investigate_target)
  ref <- ref_data[[64]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[65,], investigate_target)
  ref <- ref_data[[65]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[66,], investigate_target)
  ref <- ref_data[[66]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[67,], investigate_target)
  ref <- ref_data[[67]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[68,], investigate_target)
  ref <- ref_data[[68]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[69,], investigate_target)
  ref <- ref_data[[69]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[70,], investigate_target)
  ref <- ref_data[[70]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[71,], investigate_target)
  ref <- ref_data[[71]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[72,], investigate_target)
  ref <- ref_data[[72]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing OE-KD Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[73,], investigate_target)
  ref <- ref_data[[73]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[74,], investigate_target)
  ref <- ref_data[[74]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[75,], investigate_target)
  ref <- ref_data[[75]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[76,], investigate_target)
  ref <- ref_data[[76]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[77,], investigate_target)
  ref <- ref_data[[77]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[78,], investigate_target)
  ref <- ref_data[[78]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[79,], investigate_target)
  ref <- ref_data[[79]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[80,], investigate_target)
  ref <- ref_data[[80]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[81,], investigate_target)
  ref <- ref_data[[81]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[82,], investigate_target)
  ref <- ref_data[[82]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[83,], investigate_target)
  ref <- ref_data[[83]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[84,], investigate_target)
  ref <- ref_data[[84]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing OE-OE Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[85,], investigate_target)
  ref <- ref_data[[85]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[86,], investigate_target)
  ref <- ref_data[[86]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[87,], investigate_target)
  ref <- ref_data[[87]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[88,], investigate_target)
  ref <- ref_data[[88]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[89,], investigate_target)
  ref <- ref_data[[89]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[90,], investigate_target)
  ref <- ref_data[[90]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[91,], investigate_target)
  ref <- ref_data[[91]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[92,], investigate_target)
  ref <- ref_data[[92]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[93,], investigate_target)
  ref <- ref_data[[93]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[94,], investigate_target)
  ref <- ref_data[[94]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[95,], investigate_target)
  ref <- ref_data[[95]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[96,], investigate_target)
  ref <- ref_data[[96]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Testing OE-CP Results

# Concordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[97,], investigate_target)
  ref <- ref_data[[97]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[98,], investigate_target)
  ref <- ref_data[[98]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[99,], investigate_target)
  ref <- ref_data[[99]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[100,], investigate_target)
  ref <- ref_data[[100]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[101,], investigate_target)
  ref <- ref_data[[101]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[102,], investigate_target)
  ref <- ref_data[[102]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity > 0))
  expect_true(all(res$Similarity > 0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

# Disoncordant Results

test_that("0.2 threshold works with no cell line", {
  res <- pmap_dfr(s[103,], investigate_target)
  ref <- ref_data[[103]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.2 threshold works with one cell line", {
  res <- pmap_dfr(s[104,], investigate_target)
  ref <- ref_data[[104]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.2 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[105,], investigate_target)
  ref <- ref_data[[105]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.2))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

test_that("0.4 threshold works with no cell line", {
  res <- pmap_dfr(s[106,], investigate_target)
  ref <- ref_data[[106]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) > 2)
})

test_that("0.4 threshold works with one cell line", {
  res <- pmap_dfr(s[107,], investigate_target)
  ref <- ref_data[[107]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 1)
})

test_that("0.4 threshold works with multiple cell lines", {
  res <- pmap_dfr(s[108,], investigate_target)
  ref <- ref_data[[108]]
  expect_equal(nrow(res), nrow(ref))
  expect_true(all(res$Similarity < 0))
  expect_true(all(res$Similarity < -0.4))
  expect_true(length(unique(res$TargetCellLine)) == 2)
})

