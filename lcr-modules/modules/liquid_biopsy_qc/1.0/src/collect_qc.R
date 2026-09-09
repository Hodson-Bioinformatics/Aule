#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(); i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i]); i <- i + 1L
    if (i > length(args)) stop("Missing value for --", key, call. = FALSE)
    out[[key]] <- args[i]; i <- i + 1L
  }
  out
}

read_tsv <- function(path) {
  read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
             comment.char = "", quote = "", stringsAsFactors = FALSE)
}

number <- function(x) suppressWarnings(as.numeric(x))

first_number <- function(x) {
  x <- number(x)
  hit <- x[is.finite(x)]
  if (length(hit)) hit[1L] else NA_real_
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(as.character(x))) y else x
}

normalise_assay_mode <- function(x) {
  mode <- tolower(as.character(x %||% "capture"))
  if (!mode %in% c("capture", "swgs"))
    stop("assay-mode must be capture or swgs", call. = FALSE)
  mode
}

kv_string <- function(value) {
  if (is.null(value) || !nzchar(value)) return(list())
  pieces <- strsplit(value, ";", fixed = TRUE)[[1L]]
  setNames(lapply(pieces, function(piece) {
    fields <- strsplit(piece, "=", fixed = TRUE)[[1L]]
    paste(fields[-1L], collapse = "=")
  }), vapply(pieces, function(piece) strsplit(piece, "=", fixed = TRUE)[[1L]][1L], ""))
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  unname(cor(x[keep], y[keep], method = "spearman"))
}

residual_lowess_abs_max <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L || length(unique(x[keep])) < 3L) return(NA_real_)
  fit <- lowess(x[keep], y[keep], f = 2 / 3, iter = 3)
  max(abs(fit$y), na.rm = TRUE)
}

at_grid_edge <- function(normal_estimate, normal_grid_max, tolerance = 0.01) {
  if (!is.finite(normal_estimate) || !is.finite(normal_grid_max)) return(NA)
  normal_estimate >= normal_grid_max - tolerance
}

add_reason <- function(reasons, metric, reason) {
  reasons[[metric]] <- reason
  reasons
}

samtools_stats <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sn <- lines[startsWith(lines, "SN\t")]
  fields <- strsplit(sn, "\t", fixed = TRUE)
  keys <- sub(":$", "", vapply(fields, function(x) trimws(x[2L]), ""))
  values <- vapply(fields, function(x) first_number(x[3L]), 0)
  setNames(as.list(values), keys)
}

picard_table <- function(path, required_column) {
  lines <- readLines(path, warn = FALSE)
  header <- which(vapply(strsplit(lines, "\t", fixed = TRUE),
                         function(x) required_column %in% x, logical(1)))[1L]
  if (is.na(header) || header == length(lines)) return(NULL)
  data <- tryCatch(read.table(
    text = paste(lines[header:length(lines)], collapse = "\n"),
    sep = "\t", header = TRUE, check.names = FALSE,
    comment.char = "#", fill = TRUE, stringsAsFactors = FALSE
  ), error = function(e) NULL)
  if (is.null(data) || !required_column %in% names(data)) NULL else data
}

insert_histogram <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- grep("^insert_size\\t", lines, ignore.case = TRUE)[1L]
  if (is.na(header)) return(NULL)
  rows <- lines[(header + 1L):length(lines)]
  rows <- rows[grepl("^[0-9]+\\t", rows)]
  if (!length(rows)) return(NULL)
  table <- read.table(text = paste(c(lines[header], rows), collapse = "\n"),
                      sep = "\t", header = TRUE, check.names = FALSE)
  data.frame(size = number(table[[1L]]), count = number(table[[2L]]))
}

weighted_fraction <- function(hist, predicate) {
  if (is.null(hist) || !sum(hist$count, na.rm = TRUE)) return(NA_real_)
  sum(hist$count[predicate(hist$size)], na.rm = TRUE) /
    sum(hist$count, na.rm = TRUE)
}

canonical_chr <- function(x) {
  toupper(sub("^CHR", "", as.character(x), ignore.case = TRUE))
}

read_wig <- function(path) {
  lines <- readLines(path, warn = FALSE)
  rows <- list(); chromosome <- NA_character_; start <- step <- span <- NA_real_
  offset <- 0L; variable <- FALSE
  field <- function(line, name) {
    hit <- regmatches(line, regexpr(paste0(name, "=[^ ]+"), line))
    if (!length(hit) || !nzchar(hit)) return(NA_character_)
    sub("^[^=]+=", "", hit)
  }
  for (line in lines) {
    if (grepl("^fixedStep", line)) {
      chromosome <- field(line, "chrom"); start <- number(field(line, "start"))
      step <- number(field(line, "step")); span <- number(field(line, "span"))
      if (!is.finite(span)) span <- step
      offset <- 0L; variable <- FALSE; next
    }
    if (grepl("^variableStep", line)) {
      chromosome <- field(line, "chrom"); span <- number(field(line, "span"))
      if (!is.finite(span)) span <- 1
      variable <- TRUE; next
    }
    if (grepl("^(track|#)", line) || !nzchar(trimws(line))) next
    pieces <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    if (variable && length(pieces) >= 2L) {
      position <- number(pieces[1L]); value <- number(pieces[2L])
    } else {
      position <- start + offset * step
      value <- number(pieces[1L]); offset <- offset + 1L
    }
    if (is.finite(position) && is.finite(value)) {
      rows[[length(rows) + 1L]] <- data.frame(
        chr = canonical_chr(chromosome), start = position,
        end = position + span - 1, value = value,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame(
    chr = character(), start = numeric(), end = numeric(), value = numeric()
  )
}

detect_col <- function(table, choices) {
  hit <- which(tolower(names(table)) %in% tolower(choices))
  if (length(hit)) names(table)[hit[1L]] else NA_character_
}

bin_frame <- function(object) {
  candidates <- list(object)
  if (is.list(object)) candidates <- c(candidates, unname(object))
  for (candidate in candidates) {
    table <- tryCatch(as.data.frame(candidate), error = function(e) NULL)
    if (is.null(table)) next
    chr <- detect_col(table, c("chr", "chrom", "chromosome", "seqnames"))
    start <- detect_col(table, c("start", "pos", "begin"))
    end <- detect_col(table, c("end", "stop"))
    if (!anyNA(c(chr, start, end))) {
      value <- detect_col(table, c("median", "value", "pon"))
      return(data.frame(
        chr = canonical_chr(table[[chr]]),
        start = number(table[[start]]), end = number(table[[end]]),
        value = if (is.na(value)) NA_real_ else number(table[[value]])
      ))
    }
  }
  NULL
}

read_corrected <- function(path) {
  table <- read_tsv(path)
  chr <- detect_col(table, c("chr", "chrom", "chromosome", "seqnames"))
  start <- detect_col(table, c("start", "pos", "begin"))
  end <- detect_col(table, c("end", "stop"))
  logr <- grep("log2_TNratio_corrected|logR|log2ratio",
               names(table), value = TRUE)[1L]
  if (anyNA(c(chr, start, end, logr)))
    stop("Corrected-depth columns not recognised", call. = FALSE)
  data.frame(chr = canonical_chr(table[[chr]]),
             start = number(table[[start]]), end = number(table[[end]]),
             logr = number(table[[logr]]))
}

parse_params <- function(path) {
  lines <- readLines(path, warn = FALSE)
  labelled <- function(label) {
    hit <- grep(paste0("^", label, ":"), lines, value = TRUE)[1L]
    if (is.na(hit)) NA_character_ else trimws(sub("^[^:]+:", "", hit))
  }
  header <- grep("^init\\b", lines)[1L]
  solutions <- NULL
  if (!is.na(header) && header < length(lines)) {
    rows <- lines[(header + 1L):length(lines)]
    rows <- rows[nzchar(trimws(rows))]
    solutions <- tryCatch(read.table(
      text = paste(c(lines[header], rows), collapse = "\n"),
      header = TRUE, check.names = FALSE, stringsAsFactors = FALSE
    ), error = function(e) NULL)
  }
  tf <- first_number(labelled("Tumor Fraction"))
  ploidy <- first_number(labelled("Ploidy"))
  frac_sub <- first_number(labelled("Fraction Genome Subclonal"))
  normal_est <- if (is.finite(tf)) 1 - tf else NA_real_
  delta_bic <- gap_per_bin <- convergence <- NA_real_
  if (!is.null(solutions)) {
    n_col <- grep("^n_est$", names(solutions), value = TRUE)[1L]
    bic_col <- grep("^BIC$", names(solutions), ignore.case = TRUE,
                    value = TRUE)[1L]
    ll_col <- grep("loglik", names(solutions), ignore.case = TRUE,
                   value = TRUE)[1L]
    if (!is.na(n_col)) {
      n_values <- number(solutions[[n_col]])
      if (!is.na(ll_col) && any(is.finite(number(solutions[[ll_col]])))) {
        winner <- which.max(number(solutions[[ll_col]]))
        normal_est <- n_values[winner]
        convergence <- mean(abs(n_values - normal_est) <= 0.05, na.rm = TRUE)
      }
    }
    if (!is.na(bic_col)) {
      bic <- sort(number(solutions[[bic_col]])[
        is.finite(number(solutions[[bic_col]]))
      ])
      if (length(bic) >= 2L) delta_bic <- bic[2L] - bic[1L]
    }
    if (!is.na(ll_col)) {
      ll <- sort(number(solutions[[ll_col]])[
        is.finite(number(solutions[[ll_col]]))
      ], decreasing = TRUE)
      if (length(ll) >= 2L) gap_per_bin <- ll[1L] - ll[2L]
    }
  }
  list(tf = tf, ploidy = ploidy, frac_sub = frac_sub,
       normal_est = normal_est, delta_bic = delta_bic,
       loglik_gap = gap_per_bin, convergence = convergence)
}

metadata_hge <- function(path, sample_id) {
  if (!nzchar(path) || !file.exists(path)) return(NA_real_)
  table <- read_tsv(path)
  if (!all(c("sample_id", "input_hge_per_ml") %in% names(table)))
    return(NA_real_)
  first_number(table$input_hge_per_ml[table$sample_id == sample_id])
}

collect_qc <- function(options) {
  sample_id <- options$sample
  skip_bam <- identical(tolower(options$`skip-bam`), "true")
  assay_mode <- normalise_assay_mode(options$`assay-mode`)
  coverage_scope <- if (assay_mode == "capture") "off_target" else "genome_wide"
  reasons <- list()
  metrics <- list(sample_id = sample_id, assay_mode = assay_mode,
                  coverage_scope = coverage_scope)

  corrected <- read_corrected(options$`corrected-depth`)
  burden <- read_tsv(options$burden)
  burden <- burden[burden$sample_id == sample_id | nrow(burden) == 1L,
                   , drop = FALSE]
  if (!nrow(burden))
    stop("Burden TSV has no row for sample_id ", sample_id, call. = FALSE)
  burden <- burden[1L, , drop = FALSE]
  params <- parse_params(options$params)
  pon <- tryCatch(bin_frame(readRDS(options$pon)), error = function(e) NULL)
  gc <- read_wig(options$`gc-wig`)
  map <- read_wig(options$`map-wig`)

  n <- nrow(corrected); finite <- is.finite(corrected$logr)
  metrics$n_bins_pon <- if (is.null(pon)) {
    NA_real_
  } else if (any(is.finite(pon$value))) {
    sum(pon$chr %in% as.character(1:22) & is.finite(pon$value))
  } else {
    nrow(pon)
  }
  metrics$n_bins_used <- sum(finite)
  metrics$n_bins_nonfinite <- sum(!finite)
  metrics$frac_bins_dropped <- if (n) sum(!finite) / n else NA_real_
  if (is.null(pon))
    reasons <- add_reason(reasons, "n_bins_pon", "pon_bin_geometry_unavailable")
  if (is.null(pon)) {
    metrics$pon_bin_match_frac <- NA_real_
    reasons <- add_reason(reasons, "pon_bin_match_frac",
                          "pon_bin_geometry_unavailable")
  } else {
    sample_keys <- paste(corrected$chr, corrected$start,
                         corrected$end, sep = ":")
    pon_keys <- paste(pon$chr, pon$start, pon$end, sep = ":")
    metrics$pon_bin_match_frac <- mean(sample_keys %in% pon_keys)
  }

  align_covariate <- function(wig, metric) {
    if (!nrow(wig)) {
      reasons <<- add_reason(reasons, metric, "wig_bin_count_mismatch")
      return(rep(NA_real_, n))
    }
    sample_keys <- paste(corrected$chr, corrected$start, sep = ":")
    wig_keys <- paste(wig$chr, wig$start, sep = ":")
    matched <- wig$value[match(sample_keys, wig_keys)]
    if (sum(is.finite(matched)) < 3L && nrow(wig) == n)
      matched <- wig$value
    if (sum(is.finite(matched)) < 3L)
      reasons <<- add_reason(reasons, metric, "wig_bin_geometry_mismatch")
    matched
  }
  gc <- align_covariate(gc, "resid_gc_cor")
  map <- align_covariate(map, "resid_map_cor")
  metrics$resid_gc_cor <- safe_cor(corrected$logr, gc)
  metrics$resid_map_cor <- safe_cor(corrected$logr, map)
  metrics$resid_gc_cor_abs_max_lowess <-
    residual_lowess_abs_max(gc, corrected$logr)

  metrics$tf_estimate <- params$tf
  metrics$ploidy_estimate <- params$ploidy
  metrics$at_edge <- at_grid_edge(
    params$normal_est, number(options$`normal-grid-max`)
  )
  metrics$loglik_gap_delta_bic <- params$delta_bic
  metrics$loglik_gap_per_bin <- params$loglik_gap / max(1, metrics$n_bins_used)
  metrics$frac_genome_subclonal <- params$frac_sub
  metrics$restart_convergence <- params$convergence

  burden_names <- c(
    "burden_amp", "sigma_noise", "sigma_noise_dlrs", "vc_fit_r2",
    "vc_curv", "S_excess", "burden_amp_bootsd", "n_seg", "seg_per_gb",
    grep("^(vc_V_w|vc_nwin_w)", names(burden), value = TRUE)
  )
  for (name in unique(burden_names)) {
    metrics[[name]] <- if (name %in% names(burden)) {
      first_number(burden[[name]])
    } else {
      NA_real_
    }
    if (!name %in% names(burden))
      reasons <- add_reason(reasons, name, "missing_burden_column")
  }
  metrics$noise_ratio <- if (
    is.finite(metrics$sigma_noise) &&
    is.finite(metrics$sigma_noise_dlrs) &&
    metrics$sigma_noise_dlrs != 0
  ) {
    metrics$sigma_noise / metrics$sigma_noise_dlrs
  } else {
    NA_real_
  }
  metrics$burden_amp_ci_width_abs <- if (
    all(c("burden_amp_lo95", "burden_amp_hi95") %in% names(burden))
  ) {
    first_number(burden$burden_amp_hi95) -
      first_number(burden$burden_amp_lo95)
  } else {
    NA_real_
  }
  metrics$logr_offset_removed <- if (
    "logr_offset_removed" %in% names(burden)
  ) {
    first_number(burden$logr_offset_removed)
  } else {
    NA_real_
  }

  bam_names <- c(
    "total_reads", "mapped_reads", "mapping_rate", "dup_rate",
    "est_library_size", "insert_median", "insert_mode", "insert_sd",
    "frac_insert_lt_150", "frac_insert_150_220", "frac_insert_gt_300",
    "nucleosomal_peak_bp", "pct_on_target", "pct_off_target",
    "fold_80_base_penalty", "passing_read_ratio", "input_hge_per_ml",
    "coverage_reads", "coverage_read_frac_observed",
    "mean_coverage_depth", "median_coverage_depth",
    "coverage_depth_cv", "orphan_read_frac",
    "mean_offtarget_depth", "median_offtarget_depth", "depth_cv",
    "offtarget_reads", "offtarget_read_frac_observed"
  )
  for (name in bam_names) metrics[[name]] <- NA_real_
  metrics$input_hge_per_ml <- metadata_hge(options$metadata, sample_id)
  if (!is.finite(metrics$input_hge_per_ml))
    reasons <- add_reason(reasons, "input_hge_per_ml",
                          "metadata_unavailable")

  if (skip_bam) {
    for (name in setdiff(bam_names, "input_hge_per_ml"))
      reasons <- add_reason(reasons, name, "skip_bam_metrics")
  } else {
    prefix <- options$`bam-prefix`
    full <- samtools_stats(paste0(prefix, ".samtools_stats.txt"))
    metrics$total_reads <- first_number(full[["raw total sequences"]])
    metrics$mapped_reads <- first_number(full[["reads mapped"]])
    metrics$mapping_rate <- metrics$mapped_reads / metrics$total_reads

    markdup <- picard_table(
      paste0(prefix, ".markduplicates.txt"), "PERCENT_DUPLICATION"
    )
    if (!is.null(markdup)) {
      metrics$dup_rate <- first_number(markdup$PERCENT_DUPLICATION)
      metrics$est_library_size <-
        first_number(markdup$ESTIMATED_LIBRARY_SIZE)
    }
    insert <- picard_table(
      paste0(prefix, ".insert_size_metrics.txt"), "MEDIAN_INSERT_SIZE"
    )
    if (!is.null(insert)) {
      metrics$insert_median <- first_number(insert$MEDIAN_INSERT_SIZE)
      metrics$insert_mode <- first_number(insert$MODE_INSERT_SIZE)
      metrics$insert_sd <- first_number(insert$STANDARD_DEVIATION)
    }
    hist <- insert_histogram(paste0(prefix, ".insert_size_metrics.txt"))
    metrics$frac_insert_lt_150 <-
      weighted_fraction(hist, function(x) x < 150)
    metrics$frac_insert_150_220 <-
      weighted_fraction(hist, function(x) x >= 150 & x <= 220)
    metrics$frac_insert_gt_300 <-
      weighted_fraction(hist, function(x) x > 300)
    if (!is.null(hist) && nrow(hist))
      metrics$nucleosomal_peak_bp <- hist$size[which.max(hist$count)]

    if (assay_mode == "capture") {
      hs <- picard_table(
        paste0(prefix, ".hs_metrics.txt"), "FOLD_80_BASE_PENALTY"
      )
      if (!is.null(hs)) {
        metrics$pct_on_target <- first_number(hs$PCT_SELECTED_BASES)
        metrics$pct_off_target <- first_number(hs$PCT_OFF_BAIT)
        metrics$fold_80_base_penalty <-
          first_number(hs$FOLD_80_BASE_PENALTY)
      }
    } else {
      for (name in c(
        "pct_on_target", "pct_off_target", "fold_80_base_penalty"
      )) {
        reasons <- add_reason(reasons, name, "not_applicable_swgs")
      }
    }
    passing <- first_number(readLines(
      paste0(prefix, ".passing_reads.txt"), warn = FALSE
    ))
    metrics$passing_read_ratio <- passing / metrics$total_reads

    coverage <- samtools_stats(options$`coverage-stats`)
    metrics$coverage_reads <-
      first_number(coverage[["raw total sequences"]])
    metrics$coverage_read_frac_observed <-
      metrics$coverage_reads / metrics$total_reads
    flagstat <- readLines(options$`coverage-flagstat`, warn = FALSE)
    proper <- first_number(sub(
      " .*", "",
      grep("properly paired", flagstat, value = TRUE)[1L]
    ))
    metrics$orphan_read_frac <- if (
      is.finite(metrics$coverage_reads) && metrics$coverage_reads > 0
    ) {
      1 - proper / metrics$coverage_reads
    } else {
      NA_real_
    }
    summary <- read_tsv(options$`coverage-summary`)
    total <- summary[summary[[1L]] == "total", , drop = FALSE]
    mean_column <- grep(
      "^mean$", names(summary), ignore.case = TRUE, value = TRUE
    )[1L]
    metrics$mean_coverage_depth <- if (is.na(mean_column)) {
      NA_real_
    } else {
      first_number(total[[mean_column]])
    }
    regions <- read.table(
      gzfile(options$`coverage-regions`), sep = "\t", header = FALSE
    )
    depths <- number(regions[[ncol(regions)]])
    finite_depths <- depths[is.finite(depths)]
    metrics$median_coverage_depth <- if (length(finite_depths)) {
      median(finite_depths)
    } else {
      NA_real_
    }
    metrics$coverage_depth_cv <- if (
      length(finite_depths) && mean(finite_depths) != 0
    ) {
      sd(finite_depths) / mean(finite_depths)
    } else {
      NA_real_
    }
    metrics$depth_cv <- metrics$coverage_depth_cv

    if (assay_mode == "capture") {
      metrics$offtarget_reads <- metrics$coverage_reads
      metrics$offtarget_read_frac_observed <-
        metrics$coverage_read_frac_observed
      metrics$mean_offtarget_depth <- metrics$mean_coverage_depth
      metrics$median_offtarget_depth <- metrics$median_coverage_depth
    } else {
      for (name in c(
        "offtarget_reads", "offtarget_read_frac_observed",
        "mean_offtarget_depth", "median_offtarget_depth"
      )) {
        reasons <- add_reason(reasons, name, "not_applicable_swgs")
      }
    }
  }

  provenance <- kv_string(options$provenance)
  metrics$pon_id <- provenance$pon_id %||% basename(options$pon)
  metrics$pon_md5 <- unname(tools::md5sum(options$pon))
  metrics$pon_n_normals <- first_number(provenance$pon_n_normals)
  metrics$gc_wig_md5 <- unname(tools::md5sum(options$`gc-wig`))
  metrics$map_wig_md5 <- unname(tools::md5sum(options$`map-wig`))
  metrics$ichorcna_commit_sha <-
    provenance$ichorcna_commit_sha %||% NA_character_
  metrics$binsize_kb <- number(options$`binsize-kb`)
  metrics$genome_build <- options$`genome-build` %||% NA_character_
  metrics$genomeStyle <- provenance$genomeStyle %||% NA_character_
  metrics$pon_panel <- provenance$pon_panel %||% NA_character_
  metrics$pon_library_prep <- provenance$pon_library_prep %||% NA_character_
  metrics$pon_cohort <- provenance$pon_cohort %||% NA_character_

  thresholds <- lapply(kv_string(options$thresholds), number)
  coverage_threshold <- if (assay_mode == "capture" &&
                             !is.null(thresholds$min_mean_offtarget_depth)) {
    thresholds$min_mean_offtarget_depth
  } else {
    thresholds$min_mean_coverage_depth
  }
  flags <- character()
  flag_if <- function(condition, name) {
    if (isTRUE(condition)) flags <<- c(flags, name)
  }
  flag_if(metrics$input_hge_per_ml < thresholds$min_input_hge_per_ml,
          "low_input")
  flag_if(metrics$dup_rate > thresholds$max_dup_rate,
          "high_duplication")
  flag_if(metrics$est_library_size < thresholds$min_est_library_size,
          "low_complexity")
  flag_if(
    metrics$insert_median < thresholds$insert_median_min_bp ||
      metrics$insert_median > thresholds$insert_median_max_bp,
    "insert_size_shifted"
  )
  flag_if(metrics$frac_insert_gt_300 > thresholds$max_frac_insert_gt_300,
          "hmw_contamination")
  flag_if(metrics$passing_read_ratio < thresholds$min_passing_read_ratio,
          "undersequenced")
  flag_if(metrics$mean_coverage_depth < coverage_threshold,
          "low_coverage_depth")
  flag_if(
    is.finite(metrics$n_bins_pon) && metrics$n_bins_used != metrics$n_bins_pon,
    "bin_set_mismatch"
  )
  flag_if(abs(metrics$resid_gc_cor) > thresholds$max_abs_resid_gc_cor,
          "residual_gc_bias")
  flag_if(abs(metrics$resid_map_cor) > thresholds$max_abs_resid_map_cor,
          "residual_map_bias")
  flag_if(metrics$pon_bin_match_frac < thresholds$min_pon_bin_match_frac,
          "pon_geometry_mismatch")
  flag_if(metrics$sigma_noise > thresholds$max_sigma_noise,
          "high_noise")
  flag_if(
    metrics$noise_ratio < thresholds$noise_ratio_min ||
      metrics$noise_ratio > thresholds$noise_ratio_max,
    "noise_ratio_outlier"
  )
  if (assay_mode == "capture") {
    tolerance <- number(options$`offtarget-fraction-tolerance`)
    flag_if(
      abs(metrics$offtarget_read_frac_observed - metrics$pct_off_target) >
        tolerance,
      "offtarget_fraction_mismatch"
    )
  }

  metrics$qc_flag <- if (length(flags)) {
    paste(unique(flags), collapse = ";")
  } else {
    "pass"
  }
  for (name in names(metrics)) {
    value <- metrics[[name]]
    if (length(value) == 1L && is.na(value) &&
        is.null(reasons[[name]])) {
      reasons <- add_reason(reasons, name, "metric_unavailable")
    }
  }
  metrics$missing_reason <- if (length(reasons)) {
    paste(paste0(names(reasons), "=", unlist(reasons)), collapse = ";")
  } else {
    "none"
  }
  as.data.frame(metrics, check.names = FALSE, stringsAsFactors = FALSE)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_args(args)
  result <- collect_qc(options)
  dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
  write.table(result, options$output, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

if (!interactive() && sys.nframe() == 0L &&
    length(commandArgs(trailingOnly = TRUE))) main()
