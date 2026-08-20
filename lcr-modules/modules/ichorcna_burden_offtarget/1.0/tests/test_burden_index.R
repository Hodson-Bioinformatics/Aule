test_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
test_dir <- dirname(normalizePath(test_file))
source(file.path(test_dir, "../src/burden_index.R"))
source(file.path(test_dir, "sim_cna.R"))

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
expect_true <- function(value, info = NULL) {
  if (length(value) != 1L || is.na(value) || !isTRUE(value)) {
    stop(if (is.null(info)) "Expected TRUE" else info, call. = FALSE)
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

temporary <- tempfile("burden-tests-")
dir.create(temporary)

test_that("column detection supports correctedDepth and prefixed cna.seg", {
  simulation <- simulate_flat(seed = 2)
  corrected <- file.path(temporary, "sample.correctedDepth.txt")
  cna <- file.path(temporary, "sample.cna.seg")
  write_corrected_depth(simulation, corrected)
  write_cna_seg(simulation, cna)
  expect_equal(read_bins(corrected)$logr_column, "log2_TNratio_corrected")
  expect_equal(read_bins(cna)$logr_column, "sample.logR")
})

test_that("column detection errors without a log ratio or with copy number only", {
  bad <- file.path(temporary, "bad.tsv")
  write.table(data.frame(chr = 1, start = 1, end = 2, value = 0), bad,
              sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(read_bins(bad), "Observed columns")
  write.table(data.frame(chr = 1, start = 1, end = 2, copy.number = 2), bad,
              sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(read_bins(bad), "log-ratio")
})

test_that("locus columns are case insensitive variants", {
  path <- file.path(temporary, "variants.tsv")
  write.table(data.frame(Chromosome = 1, Begin = 1, Stop = 2, logR = 0),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_equal(read_bins(path)$data$chr, "1")
})

test_that("preparation filters, recentres, orders and counts", {
  path <- file.path(temporary, "prepare.tsv")
  table <- data.frame(
    CHROM = c("chr2", "chr1", "chr1", "chrX", "chr1"),
    POS = c(20, 30, 10, 10, 20),
    END = c(30, 40, 20, 20, 30),
    logR = c(1.2, 1.1, Inf, 1.3, 0.9)
  )
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
  prepared <- prepare_bins(path, min_bins = 2)
  expect_equal(prepared$n_nonfinite, 1)
  expect_equal(prepared$bins$chr, c("1", "1", "2"))
  expect_equal(prepared$bins$start, c(20, 30, 20))
  expect_equal(median(prepared$bins$logr), 0)
  expect_equal(prepared$offset_removed, 1.1)
  expect_error(prepare_bins(path, min_bins = 4), "min_bins")
})

test_that("exclude BED removes overlaps", {
  simulation <- simulate_flat(seed = 3)
  path <- file.path(temporary, "exclude.cna.seg")
  bed <- file.path(temporary, "exclude.bed")
  write_cna_seg(simulation, path)
  write.table(data.frame("chr1", 0, 1e7), bed, sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  unfiltered <- prepare_bins(path)
  filtered <- prepare_bins(path, bed)
  expect_true(nrow(filtered$bins) < nrow(unfiltered$bins))
  expect_true(filtered$n_excluded > 0)
})

test_that("within chromosome differences exclude boundary jumps", {
  expect_equal(within_chrom_diffs(c("1", "1", "2", "2"), c(1, 3, 100, 104)),
               c(2, 4))
})

test_that("window means respect chromosomes and trailing bins", {
  chr <- c(rep("1", 5), rep("2", 4))
  values <- 1:9
  expect_equal(window_means(chr, values, 2), c(1.5, 3.5, 6.5, 8.5))
  expect_equal(window_means(chr, values, 1), values)
})

test_that("winsorisation clips both tails and zero is identity", {
  values <- c(-100, 1:20, 100)
  clipped <- winsorise(values, 0.1)
  expect_equal(length(clipped), length(values))
  expect_true(min(clipped) > min(values) && max(clipped) < max(values))
  expect_equal(winsorise(values, 0), values)
})

test_that("variance components recover pure-noise sigma", {
  flat <- simulate_flat(sigma = 0.10, seed = 10)
  components <- variance_components(flat$chr, flat$logR, winsor_p = 0)
  expect_true(abs(sqrt(max(0, components$sigma2)) - 0.10) < 0.01)
  expect_true(abs(components$S) < 0.002)
})

test_that("variance components recover structural amplitude", {
  signal <- simulate_cna(f = 0.30, event_count = 12, sigma = 0.08, seed = 14)
  amplitude <- burden_amplitude(signal$chr, signal$logR, winsor_p = 0)
  expect_true(abs(amplitude$amp - truth_amplitude(signal)) < 0.02)
})

test_that("exported window values reproduce signed burden variance", {
  signal <- simulate_cna(f = 0.25, event_count = 12, sigma = 0.10, seed = 16)
  path <- file.path(temporary, "vc-export.cna.seg")
  write_cna_seg(signal, path)
  result <- compute_burden_index(path, n_perm = 3, n_boot = 3)
  widths <- as.numeric(sub("vc_V_w", "", grep(
    "^vc_V_w", names(result), value = TRUE
  )))
  variances <- as.numeric(result[1, paste0("vc_V_w", widths)])
  counts <- as.numeric(result[1, paste0("vc_nwin_w", widths)])
  refit <- lm(variances ~ I(1 / widths), weights = counts)
  expect_equal(unname(coef(refit)[1L]), result$burden_var_signed,
               tolerance = 1e-12)
})

test_that("diagnostic column names round-trip custom window widths", {
  signal <- simulate_cna(f = 0.25, event_count = 12, sigma = 0.10, seed = 17)
  path <- file.path(temporary, "vc-custom.cna.seg")
  write_cna_seg(signal, path)
  widths <- c(1, 3, 9, 27)
  result <- compute_burden_index(
    path, window_widths = widths, n_perm = 3, n_boot = 3
  )
  expect_equal(result$window_widths, paste(widths, collapse = ","))
  expect_true(all(paste0("vc_V_w", widths) %in% names(result)))
  expect_true(all(paste0("vc_nwin_w", widths) %in% names(result)))
})

test_that("variance components return NA on small input", {
  components <- variance_components(rep("1", 25), rnorm(25))
  expect_true(is.na(components$S))
})

test_that("amplitude is exactly invariant to ploidy offsets", {
  signal <- simulate_cna(f = 0.2, event_count = 12, sigma = 0.1, seed = 21)
  base <- burden_amplitude(signal$chr, signal$logR)$amp
  shifted <- burden_amplitude(signal$chr, signal$logR + log2(3 / 2))$amp
  expect_equal(base, shifted, tolerance = 1e-10)
})

test_that("flat-genome amplitude remains near zero across noise levels", {
  amplitudes <- vapply(c(0.05, 0.15, 0.25), function(sigma) {
    flat <- simulate_flat(sigma = sigma, seed = 33)
    burden_amplitude(flat$chr, flat$logR)$amp
  }, numeric(1))
  expect_true(all(amplitudes < 0.04))
})

test_that("amplitude is monotone in tumour fraction", {
  fractions <- c(0.02, 0.05, 0.10, 0.20, 0.40)
  amplitudes <- vapply(fractions, function(fraction) {
    signal <- simulate_cna(fraction, 12, 0.08, seed = 40)
    burden_amplitude(signal$chr, signal$logR)$amp
  }, numeric(1))
  expect_true(all(diff(amplitudes) > 0))
})

test_that("signed variance is finite and noise estimators agree", {
  flat <- simulate_flat(sigma = 0.12, seed = 55)
  amplitude <- burden_amplitude(flat$chr, flat$logR, winsor_p = 0)
  expect_true(is.finite(amplitude$var_signed))
  expect_true(abs(amplitude$sigma_noise - amplitude$sigma_noise_dlrs) < 0.02)
})

test_that("detection distinguishes flat and strong genomes", {
  flat <- simulate_flat(sigma = 0.12, seed = 60)
  signal <- simulate_cna(0.4, 20, 0.08, seed = 61)
  flat_test <- permutation_test(flat$chr, flat$logR, n_perm = 39, seed = 3)
  signal_test <- permutation_test(signal$chr, signal$logR, n_perm = 39, seed = 3)
  expect_true(flat_test$p > 0.1)
  expect_equal(signal_test$p, 1 / 40)
})

test_that("shuffling preserves total variance but changes detection statistic", {
  signal <- simulate_cna(0.3, 12, 0.08, seed = 70)
  set.seed(1)
  shuffled <- permute_within_chromosome(signal$chr, signal$logR)
  expect_equal(var(signal$logR), var(shuffled), tolerance = 1e-12)
  expect_true(abs(detect_statistic(signal$chr, signal$logR) -
                  detect_statistic(signal$chr, shuffled)) > 0.001)
})

test_that("block bootstrap is ordered and handles excessive blocks", {
  signal <- simulate_cna(0.3, 12, 0.08, seed = 80)
  estimate <- burden_amplitude(signal$chr, signal$logR)$amp
  bootstrap <- bootstrap_amplitude(signal$chr, signal$logR, n_boot = 40,
                                   block = 20, seed = 4)
  expect_true(bootstrap$lo95 <= bootstrap$hi95)
  expect_true(bootstrap$lo95 <= estimate && bootstrap$hi95 >= estimate)
  empty <- bootstrap_amplitude(signal$chr, signal$logR, n_boot = 5,
                               block = 1000)
  expect_true(is.na(empty$lo95) && is.na(empty$hi95))
})

test_that("end-to-end output is deterministic and satisfies the contract", {
  signal <- simulate_cna(0.35, 16, 0.08, seed = 90)
  cna <- file.path(temporary, "end.cna.seg")
  corrected <- file.path(temporary, "end.correctedDepth.txt")
  seg <- file.path(temporary, "end.seg.txt")
  write_cna_seg(signal, cna, "end")
  write_corrected_depth(signal, corrected)
  write_seg_txt(signal, seg, "end")
  first <- compute_burden_index(cna, "end", seg, pon_path = "/x/PON.rds",
                                n_perm = 19, n_boot = 20, seed = 7)
  second <- compute_burden_index(cna, "end", seg, pon_path = "/x/PON.rds",
                                 n_perm = 19, n_boot = 20, seed = 7)
  expect_equal(first, second, tolerance = 0)
  first_path <- file.path(temporary, "deterministic-first.tsv")
  second_path <- file.path(temporary, "deterministic-second.tsv")
  write.table(first, first_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(second, second_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_true(identical(
    readBin(first_path, "raw", n = file.info(first_path)$size),
    readBin(second_path, "raw", n = file.info(second_path)$size)
  ))
  corrected_result <- compute_burden_index(
    corrected, "end", n_perm = 19, n_boot = 20, seed = 7
  )
  expect_equal(first$burden_amp, corrected_result$burden_amp, tolerance = 1e-12)
  contract <- c(
    "sample_id", "burden_index_version", "burden_amp",
    "burden_var_signed", "burden_amp_lo95", "burden_amp_hi95",
    "burden_amp_bootsd", "detect_p", "detect_lod95", "detectable",
    "burden_snr", "sigma_noise", "sigma_noise_dlrs", "vc_fit_r2",
    "vc_curv", "S_fine", "S_coarse", "S_excess", "S_unweighted",
    "amp_coarse",
    "sd_winsor", "sd_mad", "fga_adaptive", "fga_fixed",
    "fga_threshold_used", "n_chrom_altered", "amp_per_altered_chrom",
    "burden_amp_seg", "n_seg", "seg_per_gb", "fga_seg", "n_bins_input",
    "n_bins_used", "n_bins_nonfinite", "n_bins_excluded",
    "logr_offset_removed", "logr_column", "pon_id", "binsize_kb",
    "winsor_p", "window_widths", "detect_w", "n_perm", "n_boot",
    "block_size", "seed", "qc_flag",
    "vc_V_w1", "vc_V_w2", "vc_V_w4", "vc_V_w8", "vc_V_w16",
    "vc_nwin_w1", "vc_nwin_w2", "vc_nwin_w4", "vc_nwin_w8",
    "vc_nwin_w16", "vc_resid_w1", "vc_resid_w2", "vc_resid_w4",
    "vc_resid_w8", "vc_resid_w16"
  )
  expect_equal(names(first), contract)
})

test_that("missing or unparseable segment files return NA descriptors", {
  missing <- segment_amplitude(file.path(temporary, "missing.seg.txt"))
  expect_true(all(is.na(unlist(missing))))
  bad <- file.path(temporary, "bad.seg.txt")
  writeLines("not\ta\tsegmentation", bad)
  invalid <- segment_amplitude(bad)
  expect_true(all(is.na(unlist(invalid))))
})

test_that("end-to-end exclusion reduces used bins", {
  signal <- simulate_cna(0.2, 12, 0.08, seed = 95)
  cna <- file.path(temporary, "end-exclude.cna.seg")
  bed <- file.path(temporary, "end-exclude.bed")
  write_cna_seg(signal, cna)
  write.table(data.frame("1", 0, 2e7), bed, sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  base <- compute_burden_index(cna, n_perm = 3, n_boot = 3)
  excluded <- compute_burden_index(cna, exclude_bed = bed,
                                   n_perm = 3, n_boot = 3)
  expect_true(excluded$n_bins_used < base$n_bins_used)
})

test_that("QC flags are assay descriptors rather than flatness flags", {
  synthetic <- list(
    amp = 0.4, sigma_noise = 0.1, sigma_noise_dlrs = 0.1,
    vc_fit_r2 = 0.95, variance_negative = FALSE
  )
  interval <- list(lo95 = 0.3, hi95 = 0.45)
  expect_equal(make_qc_flags(2000, synthetic, interval, 500), "pass")
  synthetic$amp <- 0
  interval <- list(lo95 = 0, hi95 = 0)
  expect_equal(make_qc_flags(2000, synthetic, interval, 500), "pass")
})

test_that("CLI parsing supports spaced, equals and boolean forms", {
  parsed <- parse_cli_args(c("--bins", "a", "--output=b", "--include_sex_chrom"))
  expect_equal(parsed$bins, "a")
  expect_equal(parsed$output, "b")
  expect_true(parsed$include_sex_chrom)
  expect_error(main(c("--bins", "a")), "required")
})

test_that("CLI writes a readable one-row TSV", {
  path <- file.path(temporary, "cli.cna.seg")
  output <- file.path(temporary, "cli.tsv")
  write_cna_seg(simulate_flat(seed = 100), path)
  main(c("--bins", path, "--output", output, "--n_perm=3", "--n_boot", "3"))
  result <- read.table(output, sep = "\t", header = TRUE, check.names = FALSE)
  expect_equal(nrow(result), 1)
})

test_that("source guard does not execute the CLI", {
  expression <- sprintf("source(%s); cat('sourced')",
                        shQuote(file.path(test_dir, "../src/burden_index.R")))
  command <- c("-e", shQuote(expression))
  output <- system2("Rscript", command, stdout = TRUE, stderr = TRUE)
  expect_true(any(output == "sourced"))
})

if (length(failures)) {
  cat("\nFailures:\n", paste0("- ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("\nAll burden-index tests passed.\n")
