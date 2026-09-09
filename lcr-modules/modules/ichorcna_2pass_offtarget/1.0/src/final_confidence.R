#!/usr/bin/env Rscript

# ==============================================================================
# final_confidence.R — Pass 2 Confidence Tier Assignment
# ==============================================================================
# Reads Pass 1 and Pass 2 ichorCNA params files together with the
# strategy TSV from pass1_qc.R. Assigns confidence tiers and estimate flags.
#
# Called by rule: _ichorcna_2pass_offtarget_confidence
# Output: confidence.tsv (tab-separated, one row per sample)

suppressPackageStartupMessages({
  library(optparse)
})

# --- Options parsing ---------------------------------------------------------

option_list <- list(
  make_option("--pass1_params", type = "character",
              help = "Path to Pass 1 .params.txt [required]"),
  make_option("--pass2_params", type = "character",
              help = "Path to Pass 2 .params.txt [required]"),
  make_option("--strategy", type = "character",
              help = "Strategy TSV from pass1_qc.R [required]"),
  make_option("--output", type = "character",
              help = "Output path for confidence.tsv [required]"),
  make_option("--sample", type = "character", default = "sample",
              help = "Sample ID [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
for (arg in c("pass1_params", "pass2_params", "strategy", "output")) {
  if (is.null(opt[[arg]])) {
    stop(sprintf("--%s is required", arg))
  }
}

# --- Helper: Parse metrics from params.txt -----------------------------------

parse_params_metrics <- function(params_file) {
  if (!file.exists(params_file)) {
    stop(sprintf("params file not found: %s", params_file))
  }
  
  lines <- readLines(params_file)
  
  # Extract Tumor Fraction
  tf_line <- grep("^Tumor Fraction:", lines, value = TRUE)
  if (length(tf_line) == 0) stop("Could not find 'Tumor Fraction:' in params file")
  pass_tf <- as.numeric(trimws(sub("^Tumor Fraction:\\s*", "", tf_line[1])))
  
  # Extract Ploidy
  ploidy_line <- grep("^Ploidy:", lines, value = TRUE)
  if (length(ploidy_line) == 0) stop("Could not find 'Ploidy:' in params file")
  pass_ploidy <- as.numeric(trimws(sub("^Ploidy:\\s*", "", ploidy_line[1])))
  
  # Extract MAD
  mad_line <- grep("^GC-Map correction MAD:", lines, value = TRUE)
  pass_MAD <- if (length(mad_line) > 0) {
    as.numeric(trimws(sub("^GC-Map correction MAD:\\s*", "", mad_line[1])))
  } else {
    NA
  }
  
  # Extract Student's t mean range
  tmean_line <- grep("^Student's t mean:", lines, value = TRUE)
  pass_t_mean_range <- if (length(tmean_line) > 0) {
    t_means <- as.numeric(
      strsplit(trimws(sub("^Student's t mean:\\s*", "", tmean_line[1])), ",\\s*")[[1]]
    )
    max(t_means) - min(t_means)
  } else {
    NA
  }
  
  # Parse solutions table
  init_idx <- grep("^init\\b", lines)[1]
  if (is.na(init_idx)) stop("Could not find solutions table in params file")
  
  init_lines <- lines[(init_idx + 1):length(lines)]
  init_lines <- init_lines[nchar(trimws(init_lines)) > 0]
  
  solutions <- list()
  for (line in init_lines) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) < 7) next
    
    loglik_val <- tryCatch(as.numeric(parts[7]), error = function(e) NA)
    if (is.na(loglik_val)) next
    
    solutions[[length(solutions) + 1]] <- data.frame(
      n_est = as.numeric(parts[2]),
      phi_est = as.numeric(parts[3]),
      loglik = loglik_val,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(solutions) == 0) stop("No valid solutions in params file")
  
  solutions_df <- do.call(rbind, solutions)
  
  # Compute QC metrics
  solutions_df <- solutions_df[order(-solutions_df$loglik), ]
  
  winner_n_est <- solutions_df$n_est[1]
  pass_loglik_gap <- if (nrow(solutions_df) >= 2) {
    solutions_df$loglik[1] - solutions_df$loglik[2]
  } else {
    Inf
  }
  
  pass_convergence <- mean(abs(solutions_df$n_est - winner_n_est) <= 0.05)
  pass_at_edge <- (winner_n_est >= 0.98)  # Close to 0.99 threshold
  
  # Ploidy cluster detection
  solutions_df$ploidy_group <- round(solutions_df$phi_est * 2) / 2
  ploidy_groups <- split(solutions_df, solutions_df$ploidy_group)
  ploidy_max_loglik <- sapply(ploidy_groups, function(g) max(g$loglik))
  ploidy_max_loglik <- sort(ploidy_max_loglik, decreasing = TRUE)
  
  pass_multiple_ploidy <- if (length(ploidy_max_loglik) >= 2) {
    (ploidy_max_loglik[1] - ploidy_max_loglik[2]) < 50
  } else {
    FALSE
  }
  
  pass_ploidy_loglik_gap <- if (length(ploidy_max_loglik) >= 2) {
    ploidy_max_loglik[1] - ploidy_max_loglik[2]
  } else {
    Inf
  }
  
  list(
    tf = pass_tf,
    ploidy = pass_ploidy,
    MAD = pass_MAD,
    t_mean_range = pass_t_mean_range,
    loglik_gap = pass_loglik_gap,
    convergence = pass_convergence,
    at_edge = pass_at_edge,
    multiple_ploidy = pass_multiple_ploidy,
    ploidy_loglik_gap = pass_ploidy_loglik_gap
  )
}

# --- Helper: Compute convergence quality label --------------------------------

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

# --- Read strategy TSV -------------------------------------------------------

strat_lines <- readLines(opt$strategy)
strat_pairs <- strsplit(strat_lines, "\t")
strat_df <- data.frame(
  key = sapply(strat_pairs, `[[`, 1),
  value = sapply(strat_pairs, `[[`, 2),
  stringsAsFactors = FALSE
)

# Create lookup
strat <- setNames(strat_df$value, strat_df$key)

strategy <- strat["strategy"]
p1_accepted <- as.logical(strat["p1_accepted"])
p1_MAD <- as.numeric(strat["p1_MAD"])
p1_t_mean_range <- as.numeric(strat["p1_t_mean_range"])
p1_loglik_gap <- as.numeric(strat["p1_loglik_gap"])
p1_convergence <- as.numeric(strat["p1_convergence"])
p1_at_edge <- as.logical(strat["p1_at_edge"])
p1_multiple_ploidy <- as.logical(strat["p1_multiple_ploidy"])
p1_ploidy_loglik_gap <- as.numeric(strat["p1_ploidy_loglik_gap"])
pass1_tf <- as.numeric(strat["pass1_tf"])
pass1_ploidy <- as.numeric(strat["pass1_ploidy"])
mad_threshold_source <- strat["mad_threshold_source"]
mad_p90_used <- as.numeric(strat["mad_p90_used"])

# --- Parse params files ------------------------------------------------------

p1 <- parse_params_metrics(opt$pass1_params)
p2 <- parse_params_metrics(opt$pass2_params)

# --- Handle ACCEPT/UNDETECTABLE: copy Pass 1 to Pass 2 metrics ---------------

if (strategy == "ACCEPT" || strategy == "UNDETECTABLE") {
  p2$tf <- p1$tf
  p2$ploidy <- p1$ploidy
  p2$MAD <- p1$MAD
  p2$t_mean_range <- p1$t_mean_range
  p2$loglik_gap <- p1$loglik_gap
  p2$convergence <- p1$convergence
  p2$at_edge <- p1$at_edge
  p2$multiple_ploidy <- p1$multiple_ploidy
  p2$ploidy_loglik_gap <- p1$ploidy_loglik_gap
}

# --- Compute Pass 2 metrics --------------------------------------------------

p2_loglik_gap <- p2$loglik_gap
p2_convergence <- p2$convergence
p2_at_edge <- p2$at_edge
pass2_tf <- p2$tf
pass2_ploidy <- p2$ploidy
tf_abs_diff <- abs(pass1_tf - pass2_tf)
conv_quality <- compute_conv_quality(p2_convergence, p2_loglik_gap)

# --- Confidence tier assignment (apply routes in order) ----------------------

tier <- "LOW"         # Default fallback
flag <- "unknown"     # Placeholder
final_tf <- pass1_tf  # Fallback
final_ploidy <- pass1_ploidy
indeterminate_high <- FALSE

# Route 1: UNDETECTABLE strategy
if (strategy == "UNDETECTABLE") {
  # Route 1: UNDETECTABLE strategy
  tier <- "HIGH"
  flag <- "undetectable_confirmed"
  final_tf <- pass1_tf
  final_ploidy <- pass1_ploidy
} else if (strategy == "ACCEPT") {
  # Route 2: ACCEPT strategy (Pass 1 accepted)
  tier <- "HIGH"
  flag <- "pass1_accepted"
  final_tf <- pass1_tf
  final_ploidy <- pass1_ploidy
} else if (p2_at_edge && p2_convergence <= 0.3) {
  # Route 3: Pass 2 failed entirely
  tier <- "UNRESOLVABLE"
  flag <- "unresolvable_no_convergence"
  final_tf <- pass1_tf
  final_ploidy <- pass1_ploidy
} else {
  # Route 4: HIGH via Pass 2
  # Category-specific thresholds (from strategy TSV, which contains the P1 category)
  # We need to infer the category from the Pass 1 QC metrics
  p1_category <- if (p1_accepted) {
    "ACCEPT"
  } else if (p1_t_mean_range < 0.08) {
    "UNDETECTABLE"
  } else if (p1_multiple_ploidy && p1_ploidy_loglik_gap < 50) {
    "PLOIDY_AMBIGUOUS"
  } else if (p1_t_mean_range < 0.15 && p1_MAD > mad_p90_used) {
    "NOISY_WEAK"
  } else if (p1_t_mean_range < 0.15 && p1_MAD <= mad_p90_used) {
    "WEAK_SIGNAL"
  } else if (p1_t_mean_range >= 0.15 && p1_MAD > mad_p90_used) {
    "NOISY_SIGNAL"
  } else {
    "CLEAN_SIGNAL"
  }
  
  # Define thresholds per category
  gap_thresholds <- list(
    NOISY_WEAK = list(gap = 15, conv = 0.6, override = FALSE),
    WEAK_SIGNAL = list(gap = 20, conv = 0.4, override = TRUE),
    CLEAN_SIGNAL = list(gap = 30, conv = 0.4, override = TRUE),
    NOISY_SIGNAL = list(gap = 25, conv = 0.4, override = TRUE),
    PLOIDY_AMBIGUOUS = list(gap = 50, conv = 0.6, override = FALSE)
  )
  
  thresholds <- gap_thresholds[[p1_category]]
  if (is.null(thresholds)) {
    # Default fallback
    thresholds <- list(gap = 20, conv = 0.4, override = FALSE)
  }
  
  # Check for HIGH tier
  high_conditions <- (
    !p2_at_edge &&
    tf_abs_diff <= 0.05 &&
    conv_quality %in% c("STRONG", "GAP_DRIVEN") &&
    p2_loglik_gap >= thresholds$gap &&
    p2_convergence >= thresholds$conv
  )
  
  # Gap override for certain categories
  if (!high_conditions && thresholds$override) {
    high_conditions <- (
      !p2_at_edge &&
      tf_abs_diff <= 0.05 &&
      p2_loglik_gap >= 80 &&
      p2_convergence >= 0.3 &&
      conv_quality %in% c("STRONG", "GAP_DRIVEN")
    )
  }
  
  if (high_conditions) {
    tier <- "HIGH"
    flag <- "pass1_pass2_agree"
    final_tf <- pass2_tf
    final_ploidy <- pass2_ploidy
  } else if (!p2_at_edge && p2_convergence >= 0.3) {
    # Route 5: MODERATE
    tier <- "MODERATE"
    flag <- if (tf_abs_diff <= 0.05) "pass2_stable" else "pass2_diverged_from_pass1"
    final_tf <- pass2_tf
    final_ploidy <- pass2_ploidy
  } else {
    # Route 6: LOW (default)
    tier <- "LOW"
    if (p2_convergence >= 0.3 && !p2_at_edge) {
      flag <- "unresolvable_pass2_best_effort"
      final_tf <- pass2_tf
      final_ploidy <- pass2_ploidy
    } else if (p2_at_edge) {
      flag <- "pass2_at_edge"
      final_tf <- pass1_tf
      final_ploidy <- pass1_ploidy
    } else {
      flag <- "pass2_low_convergence"
      final_tf <- pass1_tf
      final_ploidy <- pass1_ploidy
    }
  }
}

# --- Post-assignment suffix flags -----------------------------------------------

# INDETERMINATE_HIGH
if (tier == "LOW" && final_tf >= 0.10 && p2_loglik_gap < 10 && 
    p2_convergence < 0.40) {
  flag <- paste0(flag, "_INDETERMINATE_HIGH")
  indeterminate_high <- TRUE
}

# CONVERGED_FLAT: demotion from HIGH to MODERATE
if (conv_quality == "TRIVIAL" && 
    strategy %in% c("CLEAN_SIGNAL", "NOISY_SIGNAL")) {
  flag <- paste0(flag, "_CONVERGED_FLAT")
  if (tier == "HIGH") {
    tier <- "MODERATE"
  }
}

# --- Write confidence.tsv ------------------------------------------------

confidence_data <- data.frame(
  sample_id = opt$sample,
  strategy = strategy,
  p1_category_accepted = as.character(p1_accepted),
  p1_MAD = round(p1_MAD, 5),
  p1_t_mean_range = round(p1_t_mean_range, 4),
  p1_loglik_gap = round(p1_loglik_gap, 2),
  p1_convergence = round(p1_convergence, 4),
  p1_at_edge = as.character(p1_at_edge),
  p1_multiple_ploidy = as.character(p1_multiple_ploidy),
  pass1_tf = round(pass1_tf, 5),
  pass1_ploidy = round(pass1_ploidy, 4),
  p2_loglik_gap = round(p2_loglik_gap, 2),
  p2_convergence = round(p2_convergence, 4),
  p2_at_edge = as.character(p2_at_edge),
  conv_quality = conv_quality,
  pass2_tf = round(pass2_tf, 5),
  pass2_ploidy = round(pass2_ploidy, 4),
  tf_abs_diff = round(tf_abs_diff, 5),
  confidence_tier = tier,
  estimate_flag = flag,
  final_tf = round(final_tf, 5),
  final_ploidy = round(final_ploidy, 4),
  indeterminate_high = as.character(indeterminate_high),
  mad_threshold_source = mad_threshold_source,
  mad_p90_used = round(mad_p90_used, 5),
  stringsAsFactors = FALSE,
  row.names = FALSE
)

# Ensure output directory exists
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Write TSV
write.table(
  confidence_data,
  opt$output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# --- Console summary ---------------------------------------------------------

cat(sprintf("Sample: %s\n", opt$sample))
cat(sprintf("Strategy: %s\n", strategy))
cat(sprintf("Confidence tier: %s\n", tier))
cat(sprintf("Estimate flag: %s\n", flag))
cat(sprintf("Final TF: %.4f (ploidy: %.3f)\n", final_tf, final_ploidy))
cat(sprintf("Pass 1: TF=%.4f, conv=%.3f, at_edge=%s, MAD=%.5f\n",
            pass1_tf, p1_convergence, p1_at_edge, p1_MAD))
cat(sprintf("Pass 2: TF=%.4f, conv=%.3f, at_edge=%s, conv_quality=%s\n",
            pass2_tf, p2_convergence, p2_at_edge, conv_quality))
