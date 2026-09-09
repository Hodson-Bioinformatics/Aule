test_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
test_dir <- dirname(normalizePath(test_file))
source(file.path(test_dir, "../src/plot_vc_diagnostic.R"))

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
expect_error <- function(code, pattern = NULL) {
  error <- tryCatch({ force(code); NULL }, error = identity)
  if (is.null(error)) stop("Expected an error", call. = FALSE)
  if (!is.null(pattern) && !grepl(pattern, error$message, ignore.case = TRUE)) {
    stop("Error did not match '", pattern, "': ", error$message, call. = FALSE)
  }
}

temporary <- tempfile("vc-plot-tests-")
dir.create(temporary)
table_path <- file.path(temporary, "burden.tsv")
burden <- data.frame(
  sample_id = c("001", "sample-b"),
  window_widths = c("1,3,9", "1,3,9"),
  burden_var_signed = c(0.012, 0.020),
  burden_amp = c(0.1095, 0.1414),
  sigma_noise = c(0.08, 0.10),
  vc_fit_r2 = c(0.97, 0.91),
  vc_V_w1 = c(0.0184, 0.0300),
  vc_V_w3 = c(0.0141, 0.0233),
  vc_V_w9 = c(0.0127, 0.0211),
  vc_nwin_w1 = c(900, 800),
  vc_nwin_w3 = c(300, 266),
  vc_nwin_w9 = c(100, 88),
  check.names = FALSE
)
write.table(burden, table_path, sep = "\t", quote = FALSE, row.names = FALSE)

test_that("custom widths reconstruct stored diagnostic vectors", {
  diagnostic <- read_vc_diagnostic(table_path, "001")
  expect_equal(diagnostic$widths, c(1, 3, 9))
  expect_equal(unname(diagnostic$variances), c(0.0184, 0.0141, 0.0127))
  expect_equal(unname(diagnostic$n_win), c(900, 300, 100))
})

test_that("stored fit values pass through without refitting", {
  diagnostic <- read_vc_diagnostic(table_path, "001")
  expect_equal(diagnostic$burden_var_signed, burden$burden_var_signed[1L])
  expect_equal(diagnostic$sigma_noise, burden$sigma_noise[1L])
  expect_equal(diagnostic$burden_amp, burden$burden_amp[1L])
  expect_equal(diagnostic$vc_fit_r2, burden$vc_fit_r2[1L])
})

test_that("old tables and unknown samples fail clearly", {
  old_path <- file.path(temporary, "old.tsv")
  old <- burden[!grepl("^vc_(V|nwin)_w", names(burden))]
  write.table(old, old_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(read_vc_diagnostic(old_path, "001"), "old burden TSV")
  expect_error(read_vc_diagnostic(table_path, "missing"), "not found")
})

test_that("one call writes one requested image", {
  output <- file.path(temporary, "diagnostic.pdf")
  diagnostic <- plot_vc_diagnostic(table_path, "sample-b", output)
  expect_true(file.exists(output) && file.info(output)$size > 0)
  expect_equal(diagnostic$sample_id, "sample-b")
  expect_equal(length(list.files(temporary, pattern = "[.]pdf$")), 1)
  expect_error(
    plot_vc_diagnostic(table_path, "sample-b", file.path(temporary, "x.svg")),
    "png or .pdf"
  )
})

test_that("png extension selects the PNG device", {
  output <- file.path(temporary, "diagnostic.png")
  plot_vc_diagnostic(table_path, "001", output)
  expect_true(file.exists(output) && file.info(output)$size > 0)
})

test_that("CLI requires table, sample and output", {
  parsed <- parse_plot_args(c("--table", "a", "--sample=b", "--out", "c.pdf"))
  expect_equal(parsed$table, "a")
  expect_equal(parsed$sample, "b")
  expect_error(main(c("--table", "a", "--sample", "b")), "required")
})

test_that("workflow publishes diagnostics and supports three input modes", {
  workflow <- paste(readLines(
    file.path(test_dir, "../ichorcna_burden_offtarget.smk")
  ), collapse = "\n")
  config <- paste(readLines(
    file.path(test_dir, "../config/default.yaml")
  ), collapse = "\n")
  expect_true(grepl(
    "rule _ichorcna_burden_offtarget_plot_vc_diagnostic:",
    workflow, fixed = TRUE
  ))
  expect_true(grepl(
    "output.output.burden_diagnostic", workflow, fixed = TRUE
  ))
  expect_true(grepl(
    "plot_vc_diagnostic: \"{MODSDIR}/src/plot_vc_diagnostic.R\"",
    config, fixed = TRUE
  ))
  expect_true(grepl("ichorcna_output_dir: \"\"", config, fixed = TRUE))
  expect_true(grepl(
    "def _get_burden_ichor_output(wildcards, output_name):",
    workflow, fixed = TRUE
  ))
  expect_true(grepl("bins = _get_burden_corr_depth", workflow, fixed = TRUE))
  expect_true(grepl("generated = getattr(", workflow, fixed = TRUE))
})

test_that("source guard does not execute the CLI", {
  expression <- sprintf("source(%s); cat('sourced')",
                        shQuote(file.path(test_dir, "../src/plot_vc_diagnostic.R")))
  command <- c("-e", shQuote(expression))
  output <- system2("Rscript", command, stdout = TRUE, stderr = TRUE)
  expect_true(any(output == "sourced"))
})

if (length(failures)) {
  cat("\nFailures:\n", paste0("- ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("\nAll variance-components plot tests passed.\n")
