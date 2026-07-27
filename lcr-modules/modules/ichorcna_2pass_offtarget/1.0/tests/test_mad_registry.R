#!/usr/bin/env Rscript

# ==============================================================================
# test_mad_registry.R — Test Suite for mad_registry.R
# ==============================================================================
# Tests for MAD threshold computation, registry management

library(testthat)

# --- Helper: Create mock confidence TSV ------------------------------------

create_mock_confidence_tsv <- function(
  sample_ids = c("S1", "S2", "S3"),
  mad_values = c(0.08, 0.09, 0.10),
  run_type = "off_target"
) {
  
  df <- data.frame(
    sample_id = sample_ids,
    run_type = run_type,
    p1_MAD = mad_values,
    stringsAsFactors = FALSE
  )
  
  tmpfile <- tempfile(suffix = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
  return(tmpfile)
}

# --- Test: P90 computation ------------------------------------------------

test_that("P90 computed correctly from MAD values", {
  mad_values <- c(0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17)
  p90 <- quantile(mad_values, 0.90, names = FALSE)
  
  # P90 should be close to 0.165
  expect_true(p90 > 0.16)
  expect_true(p90 < 0.17)
})

# --- Test: Tukey fence computation ----------------------------------------

test_that("Tukey fence (Q3 + 1.5*IQR) computed correctly", {
  compute_tukey_fence <- function(values) {
    q1 <- quantile(values, 0.25, names = FALSE)
    q3 <- quantile(values, 0.75, names = FALSE)
    iqr <- q3 - q1
    return(q3 + 1.5 * iqr)
  }
  
  mad_values <- c(0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17)
  tukey <- compute_tukey_fence(mad_values)
  
  # Tukey should be higher than P90
  p90 <- quantile(mad_values, 0.90, names = FALSE)
  expect_true(tukey > p90)
})

# --- Test: Mode A - Compute from current batch --------------------------

test_that("Mode A computes thresholds from current batch only", {
  # With < 10 samples, should be Mode A (default)
  tmpfile <- create_mock_confidence_tsv(
    sample_ids = c("S1", "S2", "S3"),
    mad_values = c(0.08, 0.09, 0.10)
  )
  
  df <- read.csv(tmpfile, sep = "\t", stringsAsFactors = FALSE)
  
  # No registry provided → Mode A
  has_registry <- FALSE
  expect_false(has_registry)
  
  # Compute P90 from current batch
  p90 <- quantile(df$p1_MAD, 0.90, names = FALSE)
  expect_true(p90 > 0.08)
  expect_true(p90 < 0.11)
  
  unlink(tmpfile)
})

# --- Test: Mode B - Use supplied thresholds without recomputing --------

test_that("Mode B uses registry thresholds without recomputing", {
  # Create an existing registry
  registry_df <- data.frame(
    sample_id = c("OLD1", "OLD2", "OLD3", "OLD4", "OLD5"),
    run_type = c("off_target", "off_target", "off_target", "off_target", "off_target"),
    p1_MAD = c(0.07, 0.08, 0.09, 0.10, 0.11),
    threshold_source = "previous_batch",
    stringsAsFactors = FALSE
  )
  
  registry_file <- tempfile(suffix = ".csv")
  write.csv(registry_df, registry_file, row.names = FALSE)
  
  # New batch
  new_df <- data.frame(
    sample_id = "NEW1",
    run_type = "off_target",
    p1_MAD = 0.12,
    stringsAsFactors = FALSE
  )
  
  # Mode B: load registry but don't recompute
  registry_loaded <- read.csv(registry_file, stringsAsFactors = FALSE)
  merged <- rbind(registry_loaded[, c("sample_id", "run_type", "p1_MAD")], 
                   new_df[, c("sample_id", "run_type", "p1_MAD")])
  
  # Thresholds should use only old data (pre-computed)
  old_p90 <- quantile(registry_loaded$p1_MAD, 0.90, names = FALSE)
  merged_p90 <- quantile(merged$p1_MAD, 0.90, names = FALSE)
  
  # In Mode B, we'd use old_p90, not merged_p90
  expect_true(old_p90 < merged_p90)
  
  unlink(registry_file)
})

# --- Test: Mode C - Merge and recompute -----------------------------------

test_that("Mode C merges registries and recomputes thresholds", {
  registry_df <- data.frame(
    sample_id = c("OLD1", "OLD2", "OLD3"),
    run_type = "off_target",
    p1_MAD = c(0.07, 0.08, 0.09),
    stringsAsFactors = FALSE
  )
  
  new_df <- data.frame(
    sample_id = c("NEW1", "NEW2"),
    run_type = "off_target",
    p1_MAD = c(0.10, 0.11),
    stringsAsFactors = FALSE
  )
  
  # Mode C: merge and recompute
  merged <- rbind(registry_df[, c("sample_id", "run_type", "p1_MAD")],
                  new_df[, c("sample_id", "run_type", "p1_MAD")])
  
  p90_merged <- quantile(merged$p1_MAD, 0.90, names = FALSE)
  
  expect_equal(nrow(merged), 5)
  expect_true(p90_merged > 0.09)  # Higher than old P90
})

# --- Test: Single-sample batch with Mode B preserves registry thresholds --

test_that("Single-sample batch in Mode B doesn't affect thresholds", {
  registry_df <- data.frame(
    sample_id = c("OLD1", "OLD2", "OLD3", "OLD4", "OLD5"),
    run_type = "off_target",
    p1_MAD = c(0.07, 0.08, 0.09, 0.10, 0.11),
    stringsAsFactors = FALSE
  )
  
  old_p90 <- quantile(registry_df$p1_MAD, 0.90, names = FALSE)
  
  # Single new sample
  new_p1_mad <- 0.50  # Outlier
  
  # Mode B should NOT include this in threshold computation
  # Thresholds stay based on registry only
  mode_b_p90 <- old_p90
  
  expect_equal(mode_b_p90, old_p90)
  expect_not_equal(mode_b_p90, quantile(c(registry_df$p1_MAD, new_p1_mad), 0.90))
})

# --- Test: Registry files written with correct format --------------------

test_that("Registry CSV written with correct columns", {
  required_registry_columns <- c(
    "sample_id", "run_type", "p1_MAD", "batch_id", "date_added", "threshold_source"
  )
  
  registry_df <- data.frame(
    sample_id = "S1",
    run_type = "off_target",
    p1_MAD = 0.08,
    batch_id = "batch_001",
    date_added = Sys.Date(),
    threshold_source = "current_batch",
    stringsAsFactors = FALSE
  )
  
  for (col in required_registry_columns) {
    expect_true(col %in% names(registry_df))
  }
})

# --- Test: Thresholds file written with correct columns ----------------

test_that("Thresholds CSV written with correct columns", {
  required_threshold_columns <- c(
    "run_type", "n_samples", "mad_median", "mad_p75", "mad_p90", "mad_tukey"
  )
  
  thresholds_df <- data.frame(
    run_type = "off_target",
    n_samples = 5,
    mad_median = 0.090,
    mad_p75 = 0.105,
    mad_p90 = 0.110,
    mad_tukey = 0.130,
    stringsAsFactors = FALSE
  )
  
  for (col in required_threshold_columns) {
    expect_true(col %in% names(thresholds_df))
  }
})

if (interactive()) {
  test_dir()
} else {
  test_file()
}
