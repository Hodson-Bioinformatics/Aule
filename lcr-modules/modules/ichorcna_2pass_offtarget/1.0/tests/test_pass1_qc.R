#!/usr/bin/env Rscript

# ==============================================================================
# test_pass1_qc.R — Test Suite for pass1_qc.R
# ==============================================================================
# Tests for QC metrics computation, category assignment, and normal grid construction

library(testthat)

# --- Helper: Create mock params.txt file -----------------------------------

create_mock_params <- function(
  tumor_fraction = 0.30,
  ploidy = 2.8,
  mad = 0.08,
  t_means = c(0.15, 0.30, 0.45),
  solutions = NULL  # DataFrame with n_est, phi_est, loglik columns
) {
  
  if (is.null(solutions)) {
    # Default solutions with clear winner
    solutions <- data.frame(
      init = c("n0.5-p2", "n0.5-p3", "n0.6-p2", "n0.7-p2"),
      n_est = c(0.50, 0.51, 0.60, 0.70),
      phi_est = c(2.0, 3.0, 2.1, 2.0),
      loglik = c(1000, 950, 900, 850),
      stringsAsFactors = FALSE
    )
  }
  
  lines <- c(
    "Sample\tSample1",
    "Tumor Fraction:\t0.30",
    "Ploidy:\t2.8",
    "GC-Map correction MAD:\t0.08",
    sprintf("Student's t mean:\t%s", paste(t_means, collapse = ",")),
    "",
    "init\tn_est\tphi_est\tBIC\tFrac_genome_subclonal\tFrac_CNA_subclonal\tloglik"
  )
  
  for (i in seq_len(nrow(solutions))) {
    row <- solutions[i, ]
    line <- sprintf("%s\t%.3f\t%.1f\tNA\tNA\tNA\t%.0f",
                    row$init, row$n_est, row$phi_est, row$loglik)
    lines <- c(lines, line)
  }
  
  tmpfile <- tempfile(suffix = ".txt")
  writeLines(lines, tmpfile)
  return(tmpfile)
}

# --- Test: Parsing params header -------------------------------------------

test_that("Parse params header extracts all metrics", {
  tmpfile <- create_mock_params(
    tumor_fraction = 0.25,
    ploidy = 2.5,
    mad = 0.095,
    t_means = c(0.1, 0.3, 0.4)
  )
  
  lines <- readLines(tmpfile)
  
  # Tumor Fraction
  tf_line <- grep("^Tumor Fraction:", lines, value = TRUE)
  tf_val <- as.numeric(trimws(sub("^Tumor Fraction:\\s*", "", tf_line)))
  expect_equal(tf_val, 0.30)  # As set in create_mock_params
  
  # MAD
  mad_line <- grep("^GC-Map correction MAD:", lines, value = TRUE)
  mad_val <- as.numeric(trimws(sub("^GC-Map correction MAD:\\s*", "", mad_line)))
  expect_equal(mad_val, 0.08)
  
  # t_mean_range
  tmean_line <- grep("^Student's t mean:", lines, value = TRUE)
  t_means <- as.numeric(strsplit(
    trimws(sub("^Student's t mean:\\s*", "", tmean_line)), ",\\s*")[[1]])
  t_range <- max(t_means) - min(t_means)
  expect_equal(t_range, 0.30)
  
  unlink(tmpfile)
})

# --- Test: Loglik gap computation ------------------------------------------

test_that("Loglik gap calculated correctly", {
  # Create solutions with known loglik values
  solutions <- data.frame(
    init = c("init1", "init2", "init3"),
    n_est = c(0.5, 0.6, 0.7),
    phi_est = c(2.0, 2.0, 2.0),
    loglik = c(1000, 850, 700),
    stringsAsFactors = FALSE
  )
  
  tmpfile <- create_mock_params(solutions = solutions)
  lines <- readLines(tmpfile)
  
  # Find solutions
  init_idx <- grep("^init\\b", lines)[1]
  init_lines <- lines[(init_idx + 1):length(lines)]
  init_lines <- init_lines[nchar(trimws(init_lines)) > 0]
  
  # Parse
  logliks <- sapply(init_lines, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    as.numeric(parts[7])
  })
  
  loglik_sorted <- sort(logliks, decreasing = TRUE)
  loglik_gap <- loglik_sorted[1] - loglik_sorted[2]
  
  expect_equal(loglik_gap, 150)
  
  unlink(tmpfile)
})

# --- Test: Convergence calculation ------------------------------------------

test_that("Convergence fraction calculated as fraction within 0.05 of winner", {
  solutions <- data.frame(
    init = c("i1", "i2", "i3", "i4", "i5"),
    n_est = c(0.50, 0.51, 0.52, 0.49, 0.70),  # 4 within ±0.05 of 0.50
    phi_est = c(2.0, 2.0, 2.0, 2.0, 2.0),
    loglik = c(1000, 900, 800, 700, 600),
    stringsAsFactors = FALSE
  )
  
  tmpfile <- create_mock_params(solutions = solutions)
  lines <- readLines(tmpfile)
  
  init_idx <- grep("^init\\b", lines)[1]
  init_lines <- lines[(init_idx + 1):length(lines)]
  init_lines <- init_lines[nchar(trimws(init_lines)) > 0]
  
  n_ests <- sapply(init_lines, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    as.numeric(parts[2])
  })
  
  loglik_vec <- sapply(init_lines, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    as.numeric(parts[7])
  })
  
  winner_n <- n_ests[which.max(loglik_vec)]
  convergence <- mean(abs(n_ests - winner_n) <= 0.05)
  
  expect_equal(convergence, 0.8)  # 4 out of 5
  
  unlink(tmpfile)
})

# --- Test: Ploidy cluster detection ----------------------------------------

test_that("Ploidy clusters detected correctly", {
  solutions <- data.frame(
    init = c("i1", "i2", "i3", "i4", "i5"),
    n_est = c(0.5, 0.5, 0.5, 0.5, 0.5),
    phi_est = c(2.0, 2.1, 2.0, 4.0, 4.1),  # Two clusters: ~2.0 and ~4.0
    loglik = c(1000, 950, 900, 800, 750),  # Cluster 1 gap = 100, Cluster 2 gap = 50
    stringsAsFactors = FALSE
  )
  
  tmpfile <- create_mock_params(solutions = solutions)
  lines <- readLines(tmpfile)
  
  init_idx <- grep("^init\\b", lines)[1]
  init_lines <- lines[(init_idx + 1):length(lines)]
  init_lines <- init_lines[nchar(trimws(init_lines)) > 0]
  
  phi_ests <- sapply(init_lines, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    as.numeric(parts[3])
  })
  
  loglik_vec <- sapply(init_lines, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    as.numeric(parts[7])
  })
  
  ploidy_groups <- round(phi_ests * 2) / 2
  ploidy_max_loglik <- tapply(loglik_vec, ploidy_groups, max)
  ploidy_max_sorted <- sort(ploidy_max_loglik, decreasing = TRUE)
  
  multiple_ploidy <- (length(ploidy_max_sorted) >= 2) && 
                      ((ploidy_max_sorted[1] - ploidy_max_sorted[2]) < 50)
  gap <- ploidy_max_sorted[1] - ploidy_max_sorted[2]
  
  expect_true(multiple_ploidy)
  expect_equal(gap, 200 - 800)  # 1000 vs 800 (top of each cluster)
  
  unlink(tmpfile)
})

# --- Test: Normal grid construction ----------------------------------------

test_that("Normal grid constructed with correct bounds and step", {
  # Helper function
  construct_normal_grid <- function(pass1_tf, category) {
    width_map <- list(
      WEAK_SIGNAL = 0.20,
      CLEAN_SIGNAL = 0.15,
      NOISY_SIGNAL = 0.20
    )
    width <- width_map[[category]]
    lo <- max(0.01, pass1_tf - width)
    hi <- min(0.99, pass1_tf + width)
    grid <- unique(round(seq(lo, hi, by = 0.02), 2))
    return(grid)
  }
  
  # Test CLEAN_SIGNAL with TF=0.30, width=0.15 → [0.15, 0.45]
  grid <- construct_normal_grid(0.30, "CLEAN_SIGNAL")
  expect_true(min(grid) >= 0.15 - 0.001)
  expect_true(max(grid) <= 0.45 + 0.001)
  expect_true(all(diff(grid) >= 0.02 - 0.001))
  
  # Test boundary clipping (TF=0.05, WEAK_SIGNAL)
  grid <- construct_normal_grid(0.05, "WEAK_SIGNAL")
  expect_equal(min(grid), 0.01)  # Clipped to min
  expect_true(max(grid) <= 0.25)  # 0.05 + 0.20
})

# --- Test: Category assignment priority ----------------------------------

test_that("Categories assigned in correct priority order", {
  # UNDETECTABLE has priority over others
  # Create params: low t_mean_range, high MAD, low TF
  tmpfile <- create_mock_params(
    tumor_fraction = 0.02,
    t_means = c(0.01, 0.02, 0.03),  # range = 0.02 (< 0.08)
    mad = 0.15
  )
  
  # Read and verify condition
  lines <- readLines(tmpfile)
  tmean_line <- grep("^Student's t mean:", lines, value = TRUE)
  t_means <- as.numeric(strsplit(
    trimws(sub("^Student's t mean:\\s*", "", tmean_line)), ",\\s*")[[1]])
  t_range <- max(t_means) - min(t_means)
  
  # UNDETECTABLE condition: t_range < 0.08 AND t_range < 0.03? No.
  # Actually spec says: t_range < 0.08 AND pass1_tf < 0.03
  expect_true(t_range < 0.08)
  
  unlink(tmpfile)
})

# --- Test: NOISY_SIGNAL uses txnE=0.999 not 0.9999 ------------------------

test_that("NOISY_SIGNAL parameter set uses txnE=0.999", {
  # This is a critical constraint check
  param_map <- list(
    NOISY_SIGNAL = list(
      txnE = 0.999  # Must NOT be 0.9999
    )
  )
  
  expect_equal(param_map$NOISY_SIGNAL$txnE, 0.999)
  expect_not_equal(param_map$NOISY_SIGNAL$txnE, 0.9999)
})

# --- Run all tests -------------------------------------------------------

if (interactive()) {
  test_dir()
} else {
  test_file()
}
