library(testthat)

source_candidates <- c(file.path("..", "src"), file.path("src"))
source_dir <- source_candidates[file.exists(source_candidates)][1L]
if (is.na(source_dir)) stop("Cannot locate module src directory")
source(file.path(source_dir, "collect_qc.R"))
source(file.path(source_dir, "aggregate_qc.R"))

test_that("burden correlation identifies confounded relative width only", {
  burden <- seq_len(30)
  table <- data.frame(
    burden_amp = burden,
    tf_estimate = burden + rep(c(-0.1, 0.1), 15),
    burden_amp_ci_width_rel = rev(burden),
    burden_amp_ci_width_abs = rep(c(-1, 0, 1), 10),
    burden_amp_bootsd = rep(c(0, 1, -1, 0, -1, 1), 5)
  )
  result <- burden_correlation_table(table, 0.30)
  expect_false(result$filter_eligible[result$metric == "burden_amp_ci_width_rel"])
  expect_true(result$filter_eligible[result$metric == "burden_amp_ci_width_abs"])
  expect_true(result$filter_eligible[result$metric == "burden_amp_bootsd"])
  expect_true(all(c("n_tf_estimate", "rho_tf_estimate") %in% names(result)))
  expect_false(any(grepl("pass1", names(result))))
})

test_that("at_edge uses the supplied normal grid maximum", {
  expect_true(at_grid_edge(0.945, 0.95))
  expect_false(at_grid_edge(0.945, 0.99))
})

test_that("Tier 3 residual functions recover monotone and non-monotone bias", {
  gc <- seq(0.3, 0.7, length.out = 101)
  monotone <- gc * 2
  curved <- (gc - 0.5)^2
  expect_equal(safe_cor(monotone, gc), 1)
  expect_lt(abs(safe_cor(curved, gc)), 0.05)
  expect_gt(residual_lowess_abs_max(gc, curved), 0.01)
})

test_that("assay mode accepts capture and sWGS case-insensitively", {
  expect_equal(normalise_assay_mode("capture"), "capture")
  expect_equal(normalise_assay_mode("sWGS"), "swgs")
})

test_that("unknown assay mode is rejected", {
  expect_error(normalise_assay_mode("targeted"), "capture or swgs")
})
