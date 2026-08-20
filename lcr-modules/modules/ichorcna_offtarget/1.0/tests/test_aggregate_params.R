test_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
test_dir <- dirname(normalizePath(test_file))
source(file.path(test_dir, "../src/aggregate_params.R"))

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

write_params <- function(path, sample_id, tumor_fraction, ploidy) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "Sample\tTumor Fraction\tPloidy",
    paste(sample_id, tumor_fraction, ploidy, sep = "\t"),
    "",
    sample_id,
    "Gender:\tunknown",
    paste0("Tumor Fraction:\t", tumor_fraction),
    paste0("Ploidy:\t", ploidy),
    "Subclone Fraction:\tNA",
    "Fraction Genome Subclonal:\t0.1",
    "Fraction CNA Subclonal:\t0.2",
    "Coverage:\tNA",
    "ChrY coverage fraction:\t0.001",
    "GC-Map correction MAD:\t0.06336"
  ), path)
}

temporary <- tempfile("ichorcna-param-tests-")
first_path <- file.path(
  temporary, "capture--hg38", "param", "bin1000kb", "sample-b.param.txt"
)
second_path <- file.path(
  temporary, "capture--hg38", "param", "bin1000kb", "sample-a.param.txt"
)
write_params(first_path, "sample-b", 0.12, 2)
write_params(second_path, "sample-a", 0.05, 2)

test_that("selected solution and labelled fields are parsed", {
  result <- read_ichorcna_params(first_path)
  expect_equal(result$sample_id, "sample-b")
  expect_equal(result$tumor_fraction, 0.12)
  expect_equal(result$ploidy, 2)
  expect_true(is.na(result$subclone_fraction))
  expect_equal(result$fraction_genome_subclonal, 0.1)
  expect_equal(result$fraction_cna_subclonal, 0.2)
  expect_equal(result$chry_coverage_fraction, 0.001)
  expect_equal(result$gc_map_correction_mad, 0.06336)
})

test_that("path metadata is retained", {
  result <- read_ichorcna_params(first_path)
  expect_equal(result$seq_type, "capture")
  expect_equal(result$genome_build, "hg38")
  expect_equal(result$binsize_kb, 1000)
})

test_that("aggregation writes one sorted cohort table", {
  output <- file.path(temporary, "aggregate", "all_samples.ichorcna_summary.tsv")
  result <- aggregate_params(temporary, output)
  expect_equal(result$sample_id, c("sample-a", "sample-b"))
  expect_true(file.exists(output))
  written <- read.table(
    output, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_equal(written$sample_id, result$sample_id)
})

test_that("malformed and duplicate params fail clearly", {
  malformed <- file.path(temporary, "malformed.txt")
  writeLines("not an ichorCNA params file", malformed)
  expect_error(read_ichorcna_params(malformed), "selected-solution")
  duplicate_dir <- file.path(temporary, "duplicates")
  write_params(file.path(duplicate_dir, "one.param.txt"), "duplicate", 0.1, 2)
  write_params(file.path(duplicate_dir, "two.param.txt"), "duplicate", 0.2, 2)
  expect_error(
    aggregate_params(duplicate_dir, file.path(duplicate_dir, "out.tsv")),
    "Duplicate sample_id"
  )
})

test_that("CLI validates required and unknown arguments", {
  parsed <- parse_aggregate_args(c("--input_dir", "a", "--output=b"))
  expect_equal(parsed$input_dir, "a")
  expect_equal(parsed$output, "b")
  expect_error(main(c("--input_dir", "a")), "required")
  expect_error(
    main(c("--input_dir", "a", "--output", "b", "--extra", "c")),
    "Unknown"
  )
})

test_that("source guard does not execute the CLI", {
  expression <- sprintf("source(%s); cat('sourced')",
                        shQuote(file.path(test_dir, "../src/aggregate_params.R")))
  command <- c("-e", shQuote(expression))
  output <- system2("Rscript", command, stdout = TRUE, stderr = TRUE)
  expect_true(any(output == "sourced"))
})

if (length(failures)) {
  cat("\nFailures:\n", paste0("- ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("\nAll ichorCNA parameter aggregation tests passed.\n")
