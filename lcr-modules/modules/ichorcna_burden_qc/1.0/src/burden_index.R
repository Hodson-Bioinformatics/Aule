BURDEN_INDEX_VERSION <- "1.0"

canonical_chromosome <- function(chr) {
  value <- toupper(sub("^CHR", "", trimws(as.character(chr)), ignore.case = TRUE))
  value[value == "23"] <- "X"
  value[value == "24"] <- "Y"
  value
}

detect_column <- function(names, patterns, label) {
  for (pattern in patterns) {
    hit <- grep(pattern, names, perl = TRUE)
    if (length(hit) > 0L) {
      return(names[hit[1L]])
    }
  }
  stop(
    sprintf("Could not detect %s column. Observed columns: %s",
            label, paste(names, collapse = ", ")),
    call. = FALSE
  )
}

detect_locus_column <- function(names, choices, label) {
  hit <- which(tolower(names) %in% choices)
  if (length(hit) == 0L) {
    stop(
      sprintf("Could not detect %s column. Observed columns: %s",
              label, paste(names, collapse = ", ")),
      call. = FALSE
    )
  }
  names[hit[1L]]
}

read_bed <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("Exclude BED does not exist: ", path, call. = FALSE)
  }
  bed <- tryCatch(
    read.table(path, sep = "\t", header = FALSE, comment.char = "#",
               quote = "", stringsAsFactors = FALSE, fill = TRUE),
    error = function(error) stop("Could not read exclude BED: ", error$message,
                                 call. = FALSE)
  )
  if (ncol(bed) < 3L) {
    stop("Exclude BED must contain at least three columns", call. = FALSE)
  }
  data.frame(
    chr = canonical_chromosome(bed[[1L]]),
    start = as.numeric(bed[[2L]]),
    end = as.numeric(bed[[3L]]),
    stringsAsFactors = FALSE
  )
}

overlaps_bed <- function(chr, start, end, bed) {
  excluded <- rep(FALSE, length(chr))
  if (is.null(bed) || nrow(bed) == 0L) {
    return(excluded)
  }
  for (chromosome in intersect(unique(chr), unique(bed$chr))) {
    bin_indices <- which(chr == chromosome)
    bed_indices <- which(bed$chr == chromosome)
    for (bed_index in bed_indices) {
      excluded[bin_indices] <- excluded[bin_indices] |
        (start[bin_indices] < bed$end[bed_index] &
         end[bin_indices] > bed$start[bed_index])
    }
  }
  excluded
}

read_bins <- function(path) {
  if (!file.exists(path)) {
    stop("Bins file does not exist: ", path, call. = FALSE)
  }
  table <- tryCatch(
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
               comment.char = "", quote = "", stringsAsFactors = FALSE),
    error = function(error) stop("Could not read bins file: ", error$message,
                                 call. = FALSE)
  )
  names_observed <- names(table)
  logr_column <- detect_column(
    names_observed,
    c("^log2_TNratio_corrected$", "\\.log2_TNratio_corrected$",
      "\\.logR$", "^logR$", "^log2ratio$"),
    "log-ratio"
  )
  chr_column <- detect_locus_column(
    names_observed, c("chr", "chrom", "chromosome", "seqnames"), "chromosome"
  )
  start_column <- detect_locus_column(
    names_observed, c("start", "pos", "begin"), "start"
  )
  end_column <- detect_locus_column(
    names_observed, c("end", "stop"), "end"
  )
  data <- data.frame(
    chr = canonical_chromosome(table[[chr_column]]),
    start = suppressWarnings(as.numeric(table[[start_column]])),
    end = suppressWarnings(as.numeric(table[[end_column]])),
    logr = suppressWarnings(as.numeric(table[[logr_column]])),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(data = data, logr_column = logr_column, n_input = nrow(table))
}

prepare_bins <- function(path, exclude_bed = NULL, include_sex_chrom = FALSE,
                         min_bins = 200L) {
  parsed <- read_bins(path)
  bins <- parsed$data
  finite_logr <- is.finite(bins$logr)
  n_nonfinite <- sum(!finite_logr)
  bins <- bins[finite_logr, , drop = FALSE]
  if (any(!is.finite(bins$start) | !is.finite(bins$end))) {
    stop("Non-finite locus coordinates remain after reading bins", call. = FALSE)
  }

  if (!include_sex_chrom) {
    bins <- bins[bins$chr %in% as.character(1:22), , drop = FALSE]
  }

  bed <- read_bed(exclude_bed)
  excluded <- overlaps_bed(bins$chr, bins$start, bins$end, bed)
  n_excluded <- sum(excluded)
  bins <- bins[!excluded, , drop = FALSE]
  chr_order <- match(bins$chr, c(as.character(1:22), "X", "Y"))
  chr_order[is.na(chr_order)] <- 100L
  bins <- bins[order(chr_order, bins$chr, bins$start, bins$end), , drop = FALSE]

  if (nrow(bins) < as.integer(min_bins)) {
    stop(sprintf("Only %d bins remain after filtering; min_bins is %d",
                 nrow(bins), as.integer(min_bins)), call. = FALSE)
  }
  offset <- median(bins$logr)
  bins$logr <- bins$logr - offset
  rownames(bins) <- NULL
  list(
    bins = bins,
    n_input = parsed$n_input,
    n_nonfinite = n_nonfinite,
    n_excluded = n_excluded,
    offset_removed = offset,
    logr_column = parsed$logr_column
  )
}

window_means <- function(chr, x, w) {
  stopifnot(length(chr) == length(x), length(w) == 1L, w >= 1L)
  w <- as.integer(w)
  if (w == 1L) {
    return(as.numeric(x))
  }
  output <- numeric()
  for (chromosome in unique(chr)) {
    values <- x[chr == chromosome]
    n_complete <- length(values) %/% w
    if (n_complete > 0L) {
      values <- values[seq_len(n_complete * w)]
      output <- c(output, rowMeans(matrix(values, ncol = w, byrow = TRUE)))
    }
  }
  output
}

winsorise <- function(x, p = 0.005) {
  if (length(x) == 0L || p == 0) {
    return(x)
  }
  if (!is.finite(p) || p < 0 || p >= 0.5) {
    stop("winsor_p must be in [0, 0.5)", call. = FALSE)
  }
  limits <- as.numeric(quantile(x, probs = c(p, 1 - p), names = FALSE,
                                type = 7, na.rm = TRUE))
  pmin(pmax(x, limits[1L]), limits[2L])
}

within_chrom_diffs <- function(chr, x) {
  if (length(x) < 2L) {
    return(numeric())
  }
  differences <- diff(x)
  differences[chr[-1L] == chr[-length(chr)]]
}

variance_components <- function(chr, x, ws = c(1, 2, 4, 8, 16),
                                winsor_p = 0.005, weighted = TRUE) {
  widths <- as.numeric(ws)
  variances <- rep(NA_real_, length(widths))
  counts <- integer(length(widths))
  for (index in seq_along(widths)) {
    means <- window_means(chr, x, widths[index])
    counts[index] <- length(means)
    if (counts[index] >= 20L) {
      variances[index] <- var(winsorise(means, winsor_p))
    }
  }
  usable <- is.finite(variances) & counts >= 20L
  if (sum(usable) < 3L) {
    return(list(S = NA_real_, sigma2 = NA_real_, fit_r2 = NA_real_,
                widths = widths, variances = variances, n_win = counts,
                residuals = rep(NA_real_, length(widths))))
  }
  fit_weights <- if (weighted) counts[usable] else rep(1, sum(usable))
  fit <- lm(
    variances[usable] ~ I(1 / widths[usable]),
    weights = fit_weights
  )
  coefficients <- coef(fit)
  fit_residuals <- rep(NA_real_, length(widths))
  fit_residuals[usable] <- unname(residuals(fit))
  list(
    S = unname(coefficients[1L]),
    sigma2 = unname(coefficients[2L]),
    fit_r2 = summary(fit)$r.squared,
    widths = widths,
    variances = variances,
    n_win = counts,
    residuals = fit_residuals
  )
}

burden_amplitude <- function(chr, x, ws = c(1, 2, 4, 8, 16),
                             winsor_p = 0.005) {
  components <- variance_components(chr, x, ws, winsor_p)
  coarse <- variance_components(chr, x, c(4, 8, 16, 32), winsor_p)
  unweighted <- variance_components(chr, x, ws, winsor_p, weighted = FALSE)
  curvature_usable <- is.finite(components$variances) &
    components$n_win >= 20L
  curvature <- if (sum(curvature_usable) >= 3L) {
    inv_w <- 1 / components$widths[curvature_usable]
    fit <- lm(
      components$variances[curvature_usable] ~ inv_w + I(inv_w^2),
      weights = components$n_win[curvature_usable]
    )
    unname(coef(fit)[3L])
  } else {
    NA_real_
  }
  differences <- within_chrom_diffs(chr, x)
  sigma_dlrs <- if (length(differences) > 0L) {
    mad(differences, center = median(differences), constant = 1.4826) / sqrt(2)
  } else {
    NA_real_
  }
  list(
    amp = if (is.finite(components$S)) sqrt(max(0, components$S)) else NA_real_,
    var_signed = components$S,
    S_fine = components$S,
    S_coarse = coarse$S,
    S_excess = components$S - coarse$S,
    S_unweighted = unweighted$S,
    amp_coarse = if (is.finite(coarse$S)) sqrt(max(0, coarse$S)) else NA_real_,
    vc_curv = curvature,
    vc_widths = components$widths,
    vc_variances = components$variances,
    vc_n_win = components$n_win,
    vc_residuals = components$residuals,
    sigma_noise = if (is.finite(components$sigma2)) {
      sqrt(max(0, components$sigma2))
    } else {
      NA_real_
    },
    sigma_noise_dlrs = sigma_dlrs,
    sd_winsor = sd(winsorise(x, winsor_p)),
    sd_mad = mad(x, center = median(x), constant = 1.4826),
    vc_fit_r2 = components$fit_r2,
    variance_negative = is.finite(components$S) && components$S < 0
  )
}

vc_diagnostic_columns <- function(amplitude) {
  widths <- amplitude$vc_widths
  values <- c(
    setNames(amplitude$vc_variances, paste0("vc_V_w", widths)),
    setNames(amplitude$vc_n_win, paste0("vc_nwin_w", widths)),
    setNames(amplitude$vc_residuals, paste0("vc_resid_w", widths))
  )
  as.data.frame(as.list(values), check.names = FALSE)
}

detect_statistic <- function(chr, x, w = 8L, winsor_p = 0.005) {
  means <- window_means(chr, x, w)
  if (length(means) < 2L) {
    return(NA_real_)
  }
  var(winsorise(means, winsor_p))
}

permute_within_chromosome <- function(chr, x) {
  output <- x
  for (chromosome in unique(chr)) {
    indices <- which(chr == chromosome)
    output[indices] <- sample(x[indices], length(indices), replace = FALSE)
  }
  output
}

permutation_test <- function(chr, x, n_perm = 500L, w = 8L, seed = 1L,
                             ws = c(1, 2, 4, 8, 16), winsor_p = 0.005) {
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) {
    return(list(p = NA_real_, lod95_amp = NA_real_, observed = NA_real_,
                null_stat = numeric(), null_amp = numeric()))
  }
  set.seed(as.integer(seed))
  observed <- detect_statistic(chr, x, w, winsor_p)
  null_stat <- null_amp <- rep(NA_real_, n_perm)
  for (index in seq_len(n_perm)) {
    permuted <- permute_within_chromosome(chr, x)
    null_stat[index] <- detect_statistic(chr, permuted, w, winsor_p)
    null_amp[index] <- burden_amplitude(chr, permuted, ws, winsor_p)$amp
  }
  valid_stat <- is.finite(null_stat)
  valid_amp <- is.finite(null_amp)
  p_value <- if (is.finite(observed) && any(valid_stat)) {
    (1 + sum(null_stat >= observed, na.rm = TRUE)) / (n_perm + 1)
  } else {
    NA_real_
  }
  lod95 <- if (any(valid_amp)) {
    as.numeric(quantile(null_amp[valid_amp], 0.95, names = FALSE, type = 7))
  } else {
    NA_real_
  }
  list(p = p_value, lod95_amp = lod95, observed = observed,
       null_stat = null_stat, null_amp = null_amp)
}

bootstrap_amplitude <- function(chr, x, n_boot = 300L, block = 20L, seed = 1L,
                                ws = c(1, 2, 4, 8, 16),
                                winsor_p = 0.005) {
  eligible <- unique(chr)[vapply(unique(chr), function(chromosome) {
    sum(chr == chromosome) >= block
  }, logical(1))]
  if (length(eligible) == 0L || n_boot < 1L) {
    return(list(lo95 = NA_real_, hi95 = NA_real_, sd = NA_real_,
                amplitudes = numeric()))
  }
  set.seed(as.integer(seed))
  amplitudes <- rep(NA_real_, as.integer(n_boot))
  for (iteration in seq_len(as.integer(n_boot))) {
    boot_chr <- character()
    boot_x <- numeric()
    for (chromosome in eligible) {
      values <- x[chr == chromosome]
      needed <- length(values)
      starts <- sample.int(length(values) - block + 1L,
                           ceiling(needed / block), replace = TRUE)
      sampled <- unlist(lapply(starts, function(start) {
        values[start:(start + block - 1L)]
      }), use.names = FALSE)[seq_len(needed)]
      boot_chr <- c(boot_chr, rep(chromosome, needed))
      boot_x <- c(boot_x, sampled)
    }
    boot_x <- boot_x - median(boot_x)
    amplitudes[iteration] <- burden_amplitude(
      boot_chr, boot_x, ws, winsor_p
    )$amp
  }
  valid <- amplitudes[is.finite(amplitudes)]
  if (length(valid) == 0L) {
    return(list(lo95 = NA_real_, hi95 = NA_real_, sd = NA_real_,
                amplitudes = amplitudes))
  }
  interval <- quantile(valid, c(0.025, 0.975), names = FALSE, type = 7)
  list(lo95 = interval[1L], hi95 = interval[2L], sd = sd(valid),
       amplitudes = amplitudes)
}

fraction_altered <- function(x, sigma_noise, k = 2.5, fixed_abs = 0.20) {
  threshold <- if (is.finite(sigma_noise)) max(k * sigma_noise, 0.05) else NA_real_
  list(
    adaptive = if (is.finite(threshold)) mean(abs(x) > threshold) else NA_real_,
    fixed = mean(abs(x) > fixed_abs),
    threshold = threshold
  )
}

count_altered_chromosomes <- function(chr, x, sigma_noise, k = 3) {
  if (!is.finite(sigma_noise)) {
    return(NA_integer_)
  }
  autosomes <- intersect(unique(chr), as.character(1:22))
  sum(vapply(autosomes, function(chromosome) {
    values <- x[chr == chromosome]
    threshold <- max(k * sigma_noise / sqrt(length(values)), 0.05)
    abs(median(values)) > threshold
  }, logical(1)))
}

segment_amplitude <- function(seg_path) {
  empty <- list(amp = NA_real_, n_seg = NA_integer_, seg_per_gb = NA_real_,
                fga = NA_real_)
  if (is.null(seg_path) || is.na(seg_path) || !nzchar(seg_path) ||
      !file.exists(seg_path)) {
    return(empty)
  }
  tryCatch({
    segments <- read.table(seg_path, sep = "\t", header = TRUE,
                           check.names = FALSE, comment.char = "",
                           quote = "", stringsAsFactors = FALSE)
    names_observed <- names(segments)
    start_column <- detect_locus_column(
      names_observed, c("start", "pos", "begin"), "segment start"
    )
    end_column <- detect_locus_column(
      names_observed, c("end", "stop"), "segment end"
    )
    value_column <- detect_column(
      names_observed,
      c("^seg\\.mean$", "\\.seg\\.mean$", "^mean$", "^logR$",
        "\\.logR$", "^log2ratio$"),
      "segment log-ratio"
    )
    starts <- as.numeric(segments[[start_column]])
    ends <- as.numeric(segments[[end_column]])
    values <- as.numeric(segments[[value_column]])
    lengths <- pmax(0, ends - starts)
    valid <- is.finite(lengths) & lengths > 0 & is.finite(values)
    lengths <- lengths[valid]
    values <- values[valid]
    if (length(values) == 0L || sum(lengths) <= 0) {
      return(empty)
    }
    centre <- weighted.mean(values, lengths)
    weighted_variance <- sum(lengths * (values - centre)^2) / sum(lengths)
    list(
      amp = sqrt(weighted_variance),
      n_seg = length(values),
      seg_per_gb = length(values) / (sum(lengths) / 1e9),
      fga = sum(lengths[abs(values - centre) > 0.20]) / sum(lengths)
    )
  }, error = function(error) empty)
}

# QC flags describe assay reliability, not tumour biology. They are covariates
# or sensitivity strata and are never grounds for excluding a sample. Exclusion
# would recreate the burden-dependent selection introduced by confidence tiers.
make_qc_flags <- function(n_bins, amplitude, bootstrap, n_perm) {
  flags <- character()
  if (n_bins < 500L) flags <- c(flags, "few_bins")
  if (is.finite(amplitude$sigma_noise) && amplitude$sigma_noise > 0.25) {
    flags <- c(flags, "high_noise")
  }
  if (!is.finite(amplitude$sigma_noise)) flags <- c(flags, "noise_unestimable")
  if (isTRUE(amplitude$variance_negative)) flags <- c(flags, "variance_negative")
  if (is.finite(amplitude$vc_fit_r2) && amplitude$vc_fit_r2 < 0.8) {
    flags <- c(flags, "vc_fit_poor")
  }
  if (is.finite(bootstrap$lo95) && is.finite(bootstrap$hi95) &&
      is.finite(amplitude$amp) &&
      bootstrap$hi95 - bootstrap$lo95 > amplitude$amp) {
    flags <- c(flags, "wide_interval")
  }
  if (is.finite(amplitude$sigma_noise) &&
      is.finite(amplitude$sigma_noise_dlrs) &&
      abs(amplitude$sigma_noise - amplitude$sigma_noise_dlrs) >
      0.25 * amplitude$sigma_noise) {
    flags <- c(flags, "noise_estimators_disagree")
  }
  if (n_perm > 0L && 1 / (n_perm + 1) >= 0.05) {
    flags <- c(flags, "perm_underpowered")
  }
  if (length(flags) == 0L) "pass" else paste(flags, collapse = ";")
}

compute_burden_index <- function(
    bins_path, sample_id = NULL, seg_path = NULL, exclude_bed = NULL,
    pon_path = NULL, binsize_kb = 1000, min_bins = 200L,
    winsor_p = 0.005, window_widths = c(1, 2, 4, 8, 16),
    detect_w = 8L, n_perm = 500L, n_boot = 300L, block = 20L,
    seed = 1L, fga_k = 2.5, include_sex_chrom = FALSE) {
  prepared <- prepare_bins(bins_path, exclude_bed, include_sex_chrom, min_bins)
  bins <- prepared$bins
  amplitude <- burden_amplitude(
    bins$chr, bins$logr, window_widths, winsor_p
  )
  permutation <- permutation_test(
    bins$chr, bins$logr, n_perm, detect_w, seed, window_widths, winsor_p
  )
  bootstrap <- bootstrap_amplitude(
    bins$chr, bins$logr, n_boot, block, seed, window_widths, winsor_p
  )
  fga <- fraction_altered(bins$logr, amplitude$sigma_noise, fga_k)
  n_chrom_altered <- count_altered_chromosomes(
    bins$chr, bins$logr, amplitude$sigma_noise
  )
  segment <- segment_amplitude(seg_path)
  if (is.null(sample_id) || !nzchar(sample_id)) {
    sample_id <- sub("\\.(correctedDepth\\.txt|cna\\.seg)$", "",
                     basename(bins_path))
  }
  pon_id <- if (is.null(pon_path) || is.na(pon_path) || !nzchar(pon_path)) {
    NA_character_
  } else {
    basename(pon_path)
  }
  detectable <- is.finite(permutation$p) && permutation$p < 0.05
  output <- data.frame(
    sample_id = sample_id,
    burden_index_version = BURDEN_INDEX_VERSION,
    burden_amp = amplitude$amp,
    burden_var_signed = amplitude$var_signed,
    burden_amp_lo95 = bootstrap$lo95,
    burden_amp_hi95 = bootstrap$hi95,
    burden_amp_bootsd = bootstrap$sd,
    detect_p = permutation$p,
    detect_lod95 = permutation$lod95_amp,
    detectable = detectable,
    burden_snr = if (is.finite(amplitude$sigma_noise) &&
                     amplitude$sigma_noise > 0) {
      amplitude$amp / amplitude$sigma_noise
    } else {
      NA_real_
    },
    sigma_noise = amplitude$sigma_noise,
    sigma_noise_dlrs = amplitude$sigma_noise_dlrs,
    vc_fit_r2 = amplitude$vc_fit_r2,
    vc_curv = amplitude$vc_curv,
    S_fine = amplitude$S_fine,
    S_coarse = amplitude$S_coarse,
    S_excess = amplitude$S_excess,
    S_unweighted = amplitude$S_unweighted,
    amp_coarse = amplitude$amp_coarse,
    sd_winsor = amplitude$sd_winsor,
    sd_mad = amplitude$sd_mad,
    fga_adaptive = fga$adaptive,
    fga_fixed = fga$fixed,
    fga_threshold_used = fga$threshold,
    n_chrom_altered = n_chrom_altered,
    amp_per_altered_chrom = if (is.finite(n_chrom_altered) &&
                                n_chrom_altered > 0) {
      amplitude$amp / n_chrom_altered
    } else {
      NA_real_
    },
    burden_amp_seg = segment$amp,
    n_seg = segment$n_seg,
    seg_per_gb = segment$seg_per_gb,
    fga_seg = segment$fga,
    n_bins_input = prepared$n_input,
    n_bins_used = nrow(bins),
    n_bins_nonfinite = prepared$n_nonfinite,
    n_bins_excluded = prepared$n_excluded,
    logr_offset_removed = prepared$offset_removed,
    logr_column = prepared$logr_column,
    pon_id = pon_id,
    binsize_kb = as.numeric(binsize_kb),
    winsor_p = winsor_p,
    window_widths = paste(window_widths, collapse = ","),
    detect_w = as.integer(detect_w),
    n_perm = as.integer(n_perm),
    n_boot = as.integer(n_boot),
    block_size = as.integer(block),
    seed = as.integer(seed),
    qc_flag = make_qc_flags(nrow(bins), amplitude, bootstrap, as.integer(n_perm)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  output <- cbind(output, vc_diagnostic_columns(amplitude))
  output
}

parse_cli_args <- function(args) {
  parsed <- list()
  index <- 1L
  while (index <= length(args)) {
    token <- args[index]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
      key <- pieces[1L]
      value <- paste(pieces[-1L], collapse = "=")
    } else {
      key <- token
      if (key == "include_sex_chrom") {
        value <- TRUE
      } else {
        if (index == length(args) || startsWith(args[index + 1L], "--")) {
          stop("Missing value for --", key, call. = FALSE)
        }
        index <- index + 1L
        value <- args[index]
      }
    }
    parsed[[key]] <- value
    index <- index + 1L
  }
  parsed
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_cli_args(args)
  missing <- c("bins", "output")[!c("bins", "output") %in% names(options)]
  if (length(missing) > 0L) {
    stop("Missing required argument(s): ",
         paste(paste0("--", missing), collapse = ", "), call. = FALSE)
  }
  numeric_option <- function(name, default) {
    if (is.null(options[[name]])) default else as.numeric(options[[name]])
  }
  widths <- if (is.null(options$window_widths)) {
    c(1, 2, 4, 8, 16)
  } else {
    as.numeric(strsplit(options$window_widths, ",", fixed = TRUE)[[1L]])
  }
  result <- compute_burden_index(
    bins_path = options$bins,
    sample_id = options$sample,
    seg_path = options$seg,
    exclude_bed = options$exclude_bed,
    pon_path = options$pon,
    binsize_kb = numeric_option("binsize_kb", 1000),
    min_bins = numeric_option("min_bins", 200),
    winsor_p = numeric_option("winsor_p", 0.005),
    window_widths = widths,
    detect_w = numeric_option("detect_w", 8),
    n_perm = numeric_option("n_perm", 500),
    n_boot = numeric_option("n_boot", 300),
    block = numeric_option("block", 20),
    seed = numeric_option("seed", 1),
    fga_k = numeric_option("fga_k", 2.5),
    include_sex_chrom = isTRUE(options$include_sex_chrom)
  )
  output_directory <- dirname(options$output)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  write.table(result, options$output, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  invisible(result)
}

if (!interactive() && sys.nframe() == 0L &&
    length(commandArgs(trailingOnly = TRUE)) > 0L) main()
