test_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
test_dir <- dirname(normalizePath(test_file))
source(file.path(test_dir, "../src/aggregate_burden.R"))

failures <- character()
test_that <- function(description, code) {
  tryCatch({
    force(code)
    cat("PASS:", description, "\n")
  }, error = function(error) {
    failures <<- c(failures, paste(description, error$message, sep = ": "))
    cat("FAIL:", description, "-", error$message, "\n")
  })
}
expect_true <- function(value) {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) {
    stop("Expected TRUE", call. = FALSE)
  }
}
expect_equal <- function(object, expected, tolerance = 1e-8) {
  comparison <- all.equal(object, expected, tolerance = tolerance,
                          check.attributes = FALSE)
  if (!isTRUE(comparison)) stop(paste(comparison, collapse = "; "), call. = FALSE)
}

temporary <- tempfile("aggregate-tests-")
dir.create(temporary)
base <- data.frame(
  sample_id = c("a", "b"),
  burden_amp = c(0.1, 0.3),
  sigma_noise = c(0.2, 0.1),
  detect_p = c(0.01, 0.20),
  detectable = c(TRUE, FALSE),
  qc_flag = c("pass", "high_noise;vc_fit_poor"),
  stringsAsFactors = FALSE
)
extra <- base[1, ]
extra$sample_id <- "c"
extra$extra_column <- 1
path_one <- file.path(temporary, "one.tsv")
path_two <- file.path(temporary, "two.tsv")
write.table(base, path_one, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(extra, path_two, sep = "\t", quote = FALSE, row.names = FALSE)

test_that("heterogeneous schemas use their intersection with a warning", {
  warned <- FALSE
  table <- withCallingHandlers(
    read_burden_tables(c(path_one, path_two)),
    warning = function(warning) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(warned)
  expect_true(!"extra_column" %in% names(table))
  expect_equal(nrow(table), 3)
})

test_that("cohort annotations are bounded and BH adjusted", {
  table <- add_cohort_annotations(read_burden_tables(path_one))
  expect_true(all(table$burden_amp_rank >= 0 & table$burden_amp_rank <= 1))
  expect_true(all(table$sigma_noise_rank >= 0 & table$sigma_noise_rank <= 1))
  expect_equal(table$detect_q, p.adjust(base$detect_p, "BH"))
  expect_equal(table$detectable_fdr10, c(TRUE, FALSE))
})

test_that("QC summary reports flag counts and detection rate", {
  table <- add_cohort_annotations(read_burden_tables(path_one))
  summary <- qc_summary(table)
  high_noise <- summary[summary$summary_type == "qc_flag" &
                        summary$metric == "high_noise", ]
  detection <- summary[summary$summary_type == "detection_rate", ]
  expect_equal(high_noise$n, 1)
  expect_equal(detection$value, 0.5)
})

test_that("aggregation writes cohort and QC tables", {
  cohort_path <- file.path(temporary, "out", "cohort.tsv")
  qc_path <- file.path(temporary, "out", "qc.tsv")
  result <- suppressWarnings(aggregate_burden(
    c(path_one, path_two), cohort_path, qc_path
  ))
  expect_true(file.exists(cohort_path) && file.exists(qc_path))
  expect_equal(nrow(result$table), 3)
  expect_equal(nrow(read.table(cohort_path, header = TRUE, sep = "\t")), 3)
})

if (length(failures)) {
  cat("\nFailures:\n", paste0("- ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("\nAll aggregation tests passed.\n")
