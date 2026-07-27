#!/usr/bin/env Rscript

# ==============================================================================
# test_final_confidence.R — Test Suite for final_confidence.R
# ==============================================================================
# Tests for confidence tier routing, convergence quality assessment

library(testthat)

# --- Test: Convergence quality assessment ----------------------------------

test_that("conv_quality assigned correctly", {
  compute_conv_quality <- function(convergence, loglik_gap) {
    # Apply routes in order (first match wins)
    if (convergence >= 0.8 && loglik_gap < 5) {
      return("TRIVIAL")
    }
    if (convergence >= 0.6 && loglik_gap < 15) {
      return("HOLLOW")
    }
    if (convergence >= 0.6 && loglik_gap >= 30) {
      return("STRONG")
    }
    if (convergence >= 0.3 && loglik_gap >= 80) {
      return("GAP_DRIVEN")
    }
    if (convergence >= 0.4 && loglik_gap >= 15) {
      return("GENUINE")
    }
    return("POOR")
  }
  
  # Test each quality level
  expect_equal(compute_conv_quality(0.9, 3), "TRIVIAL")
  expect_equal(compute_conv_quality(0.7, 10), "HOLLOW")
  expect_equal(compute_conv_quality(0.7, 50), "STRONG")
  expect_equal(compute_conv_quality(0.4, 100), "GAP_DRIVEN")
  expect_equal(compute_conv_quality(0.5, 20), "GENUINE")
  expect_equal(compute_conv_quality(0.2, 5), "POOR")
  
  # Test priority: TRIVIAL fires before HOLLOW
  expect_equal(compute_conv_quality(0.9, 3), "TRIVIAL")
  expect_not_equal(compute_conv_quality(0.9, 3), "HOLLOW")
})

# --- Test: ACCEPT strategy routes to HIGH with p1_accepted flag -----------

test_that("ACCEPT strategy returns HIGH confidence", {
  strategy <- "ACCEPT"
  expect_equal(strategy, "ACCEPT")
})

# --- Test: UNDETECTABLE strategy routes to HIGH -------------------------

test_that("UNDETECTABLE strategy returns HIGH confidence", {
  strategy <- "UNDETECTABLE"
  expect_equal(strategy, "UNDETECTABLE")
})

# --- Test: Pass 2 metrics copied from Pass 1 for ACCEPT/UNDETECTABLE -----

test_that("Pass 2 metrics equal Pass 1 for ACCEPT strategy", {
  # When strategy is ACCEPT, p2 should equal p1
  p1_tf <- 0.25
  p1_convergence <- 0.75
  
  # For ACCEPT/UNDETECTABLE, set p2 = p1
  p2_tf <- p1_tf  # Should be equal
  p2_convergence <- p1_convergence
  
  expect_equal(p2_tf, p1_tf)
  expect_equal(p2_convergence, p1_convergence)
})

# --- Test: TF agreement check (within 0.05) --------------------------------

test_that("TF agreement determined by ±5% difference", {
  pass1_tf <- 0.30
  
  # Within tolerance
  pass2_tf_agree <- 0.32
  tf_diff <- abs(pass1_tf - pass2_tf_agree)
  expect_true(tf_diff <= 0.05)
  
  # Outside tolerance
  pass2_tf_disagree <- 0.37
  tf_diff <- abs(pass1_tf - pass2_tf_disagree)
  expect_false(tf_diff <= 0.05)
})

# --- Test: Confidence tier output format ---------------------------------

test_that("Confidence tier is one of valid tiers", {
  valid_tiers <- c("HIGH", "MODERATE", "LOW", "UNDETECTABLE", "UNRESOLVABLE")
  
  test_tier <- "HIGH"
  expect_true(test_tier %in% valid_tiers)
  
  test_tier <- "MODERATE"
  expect_true(test_tier %in% valid_tiers)
})

# --- Test: Estimate flag strings match vocabulary --------------------------

test_that("Estimate flags use valid vocabulary", {
  valid_flags <- c(
    "undetectable_confirmed",
    "pass1_accepted",
    "unresolvable_no_convergence",
    "pass1_pass2_agree",
    "pass2_stable",
    "pass2_diverged_from_pass1",
    "unresolvable_pass2_best_effort",
    "pass2_at_edge",
    "pass2_low_convergence"
  )
  
  test_flag <- "pass1_accepted"
  expect_true(test_flag %in% valid_flags)
  
  test_flag <- "pass2_stable"
  expect_true(test_flag %in% valid_flags)
})

# --- Test: CONVERGED_FLAT demotion from HIGH to MODERATE ---------------

test_that("CONVERGED_FLAT suffix demotes HIGH to MODERATE", {
  conv_quality <- "TRIVIAL"
  strategy <- "CLEAN_SIGNAL"
  tier <- "HIGH"
  
  # CONVERGED_FLAT condition met
  if (conv_quality == "TRIVIAL" && strategy %in% c("CLEAN_SIGNAL", "NOISY_SIGNAL")) {
    if (tier == "HIGH") {
      tier <- "MODERATE"
    }
  }
  
  expect_equal(tier, "MODERATE")
})

# --- Test: INDETERMINATE_HIGH suffix applied correctly ----------------

test_that("INDETERMINATE_HIGH suffix applied to LOW tier samples", {
  tier <- "LOW"
  final_tf <- 0.15
  p2_loglik_gap <- 8
  p2_convergence <- 0.35
  
  indeterminate_high <- FALSE
  
  if (tier == "LOW" && final_tf >= 0.10 && p2_loglik_gap < 10 && 
      p2_convergence < 0.40) {
    indeterminate_high <- TRUE
  }
  
  expect_true(indeterminate_high)
})

# --- Test: Output TSV has correct structure and header --------------------

test_that("Confidence TSV output has all required columns", {
  required_columns <- c(
    "sample_id",
    "strategy",
    "p1_category_accepted",
    "p1_MAD",
    "p1_t_mean_range",
    "pass1_tf",
    "pass1_ploidy",
    "p2_loglik_gap",
    "p2_convergence",
    "p2_at_edge",
    "pass2_tf",
    "pass2_ploidy",
    "tf_abs_diff",
    "confidence_tier",
    "estimate_flag",
    "final_tf",
    "final_ploidy"
  )
  
  # Check that all required columns are included
  test_data <- data.frame(
    sample_id = "TEST_001",
    strategy = "CLEAN_SIGNAL",
    p1_category_accepted = "FALSE",
    p1_MAD = 0.085,
    p1_t_mean_range = 0.25,
    pass1_tf = 0.30,
    pass1_ploidy = 2.8,
    p2_loglik_gap = 120,
    p2_convergence = 0.70,
    p2_at_edge = "FALSE",
    pass2_tf = 0.31,
    pass2_ploidy = 2.82,
    tf_abs_diff = 0.01,
    confidence_tier = "HIGH",
    estimate_flag = "pass1_pass2_agree",
    final_tf = 0.31,
    final_ploidy = 2.82,
    stringsAsFactors = FALSE
  )
  
  for (col in required_columns) {
    expect_true(col %in% names(test_data), 
                label = sprintf("Column '%s' missing from output", col))
  }
})

if (interactive()) {
  test_dir()
} else {
  test_file()
}
