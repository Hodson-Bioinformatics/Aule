#!/usr/bin/env Rscript

# ==============================================================================
# pass1_qc.R — Pass 1 QC and Strategy Decision Logic
# ==============================================================================
# Parses ichorCNA Pass 1 params.txt, computes QC metrics, and writes a
# strategy TSV that guides Pass 2 parameters or signals acceptance/rejection.
#
# Called by rule: _ichorcna_2pass_offtarget_pass1_qc
# Output: strategy.tsv (tab-separated key-value pairs)

suppressPackageStartupMessages({
  library(optparse)
})


# --- Options parsing ---------------------------------------------------------

option_list <- list(
  make_option("--params", type = "character",
              help = "Path to Pass 1 .params.txt [required]"),
  make_option("--strategy", type = "character",
              help = "Output path for strategy.tsv [required]"),
  make_option("--sample", type = "character", default = "sample",
              help = "Sample ID [default: %default]"),
  make_option("--max_normal_init", type = "numeric", default = 0.99,
              help = "Highest normal fraction initialized in Pass 1 [default: %default]"),
  make_option("--mad_registry_file", type = "character", default = NULL,
              help = "Path to MAD registry CSV (optional)"),
  make_option("--mad_registry_run_type", type = "character", default = "off_target",
              help = "Run type for registry lookup: 'off_target' or 'all_target' [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$params)) {
  stop("--params is required")
}
if (is.null(opt$strategy)) {
  stop("--strategy is required")
}

# --- Helper: Parse params.txt header block -----------------------------------

parse_params_header <- function(params_file) {
  if (!file.exists(params_file)) {
    stop(sprintf("params file not found: %s", params_file))
  }
  
  lines <- readLines(params_file)
  
  # Extract Tumor Fraction
  tf_line <- grep("^Tumor Fraction:", lines, value = TRUE)
  if (length(tf_line) == 0) {
    stop("Could not find 'Tumor Fraction:' in params file")
  }
  pass1_tf <- as.numeric(trimws(sub("^Tumor Fraction:\\s*", "", tf_line[1])))
  
  # Extract Ploidy
  ploidy_line <- grep("^Ploidy:", lines, value = TRUE)
  if (length(ploidy_line) == 0) {
    stop("Could not find 'Ploidy:' in params file")
  }
  pass1_ploidy <- as.numeric(trimws(sub("^Ploidy:\\s*", "", ploidy_line[1])))
  
  # Extract GC-Map correction MAD
  mad_line <- grep("^GC-Map correction MAD:", lines, value = TRUE)
  if (length(mad_line) == 0) {
    stop("Could not find 'GC-Map correction MAD:' in params file")
  }
  p1_MAD <- as.numeric(trimws(sub("^GC-Map correction MAD:\\s*", "", mad_line[1])))
  
  # Extract Student's t mean (comma-separated)
  tmean_line <- grep("^Student's t mean:", lines, value = TRUE)
  if (length(tmean_line) == 0) {
    stop("Could not find 'Student's t mean:' in params file")
  }
  t_means <- as.numeric(
    strsplit(trimws(sub("^Student's t mean:\\s*", "", tmean_line[1])), ",\\s*")[[1]]
  )
  p1_t_mean_range <- max(t_means) - min(t_means)
  
  list(
    pass1_tf = pass1_tf,
    pass1_ploidy = pass1_ploidy,
    p1_MAD = p1_MAD,
    p1_t_mean_range = p1_t_mean_range,
    lines = lines
  )
}

# --- Helper: Parse solutions table -------------------------------------------

parse_solutions_table <- function(lines) {
  # Find the "init\tn_est\tphi_est\t..." header line
  init_idx <- grep("^init\\b", lines)[1]
  if (is.na(init_idx)) {
    stop("Could not find initialization table header in params file")
  }
  
  # Extract solution rows
  init_lines <- lines[(init_idx + 1):length(lines)]
  init_lines <- init_lines[nchar(trimws(init_lines)) > 0]
  
  # Parse each row: init n_est phi_est BIC Frac_genome_subclonal Frac_CNA_subclonal loglik
  solutions <- list()
  for (line in init_lines) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) < 7) next
    
    # Protect against NA or non-numeric values
    loglik_val <- tryCatch(as.numeric(parts[7]), error = function(e) NA)
    if (is.na(loglik_val)) next
    
    solutions[[length(solutions) + 1]] <- data.frame(
      init = parts[1],
      n_est = as.numeric(parts[2]),
      phi_est = as.numeric(parts[3]),
      loglik = loglik_val,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(solutions) == 0) {
    stop("No valid solutions parsed from table")
  }
  
  do.call(rbind, solutions)
}

# --- Helper: Compute QC metrics from solutions table --------------------------

compute_qc_metrics <- function(solutions_df, max_normal_init) {
  # Sort by loglik descending
  solutions_df <- solutions_df[order(-solutions_df$loglik), ]
  
  # Extract metrics
  winner_idx <- 1
  winner_n_est <- solutions_df$n_est[winner_idx]
  winner_phi <- solutions_df$phi_est[winner_idx]
  
  # Loglik gap: difference between top two solutions
  p1_loglik_gap <- if (nrow(solutions_df) >= 2) {
    solutions_df$loglik[1] - solutions_df$loglik[2]
  } else {
    Inf
  }
  
  # Convergence: fraction of solutions with n_est within ±0.05 of winner
  p1_convergence <- mean(abs(solutions_df$n_est - winner_n_est) <= 0.05)
  
  # At edge: winner n_est close to max_normal_init
  p1_at_edge <- (winner_n_est >= (max_normal_init - 0.01))
  
  # Ploidy cluster detection: round phi_est, group, take max loglik per group
  solutions_df$ploidy_group <- round(solutions_df$phi_est * 2) / 2  # Round to nearest 0.5
  ploidy_groups <- split(solutions_df, solutions_df$ploidy_group)
  ploidy_max_loglik <- sapply(ploidy_groups, function(g) max(g$loglik))
  ploidy_max_loglik <- sort(ploidy_max_loglik, decreasing = TRUE)
  
  p1_multiple_ploidy <- if (length(ploidy_max_loglik) >= 2) {
    (ploidy_max_loglik[1] - ploidy_max_loglik[2]) < 50
  } else {
    FALSE
  }
  
  p1_ploidy_loglik_gap <- if (length(ploidy_max_loglik) >= 2) {
    ploidy_max_loglik[1] - ploidy_max_loglik[2]
  } else {
    Inf
  }
  
  list(
    winner_n_est = winner_n_est,
    winner_phi = winner_phi,
    p1_loglik_gap = p1_loglik_gap,
    p1_convergence = p1_convergence,
    p1_at_edge = p1_at_edge,
    p1_multiple_ploidy = p1_multiple_ploidy,
    p1_ploidy_loglik_gap = p1_ploidy_loglik_gap
  )
}

# --- Helper: Load MAD threshold from registry --------------------------------

get_mad_threshold <- function(registry_path, run_type) {
  # Hardcoded fallbacks
  fallback_thresholds <- list(
    off_target = 0.10,
    all_target = 0.13
  )
  
  threshold_source <- "fallback"
  mad_p90 <- fallback_thresholds[[run_type]]
  
  if (!is.null(registry_path) && file.exists(registry_path)) {
    tryCatch({
      registry <- read.csv(registry_path, stringsAsFactors = FALSE)
      
      # Filter for samples of this run_type and check count
      run_samples <- registry[registry$run_type == run_type, ]
      if (nrow(run_samples) >= 10) {
        # Load MAD P90 from registry
        mad_p90 <- quantile(run_samples$p1_MAD, 0.90, names = FALSE)
        threshold_source <- "registry"
      }
    }, error = function(e) {
      warning(sprintf("Could not load MAD registry: %s", e$message))
    })
  }
  
  if (threshold_source == "fallback") {
    warning(sprintf(
      "MAD threshold: using fallback %.3f for %s (registry has <10 samples or unavailable)",
      mad_p90, run_type
    ))
  }
  
  list(mad_p90 = mad_p90, source = threshold_source)
}

# --- Helper: Normal grid construction ----------------------------------------

construct_normal_grid <- function(pass1_tf, pass1_ploidy, category, 
                                  winner_n_est) {
  if (category == "PLOIDY_AMBIGUOUS") {
    # Broad grid
    grid <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)
  } else {
    # Category-specific width
    width_map <- list(
      UNDETECTABLE = 0.15,
      NOISY_WEAK = 0.15,
      WEAK_SIGNAL = 0.20,
      CLEAN_SIGNAL = 0.15,
      NOISY_SIGNAL = 0.20,
      PLOIDY_AMBIGUOUS = NA  # Handled above
    )
    
    width <- width_map[[category]]
    if (is.na(width)) {
      width <- 0.15  # Default
    }
    
    lo <- max(0.01, pass1_tf - width)
    hi <- min(0.99, pass1_tf + width)
    
    grid <- unique(round(seq(lo, hi, by = 0.02), 2))
  }
  
  # Format as R vector
  paste0("c(", paste(grid, collapse = ","), ")")
}

# --- Main execution ----------------------------------------------------------

# Parse params file
params <- parse_params_header(opt$params)
solutions_df <- parse_solutions_table(params$lines)

# Compute QC metrics
metrics <- compute_qc_metrics(solutions_df, opt$max_normal_init)

# Get MAD threshold
mad_info <- get_mad_threshold(opt$mad_registry_file, opt$mad_registry_run_type)
mad_p90 <- mad_info$mad_p90
mad_threshold_source <- mad_info$source

# --- Pass 1 accept check -----------------------------------------------------

p1_accepted <- (
  !metrics$p1_at_edge &&
  metrics$p1_loglik_gap >= 100 &&
  params$p1_t_mean_range >= 0.20 &&
  params$p1_MAD < mad_p90 &&
  metrics$p1_convergence >= 0.60
)

# --- Category assignment (priority order) ------------------------------------

strategy <- "CLEAN_SIGNAL"  # Default
category <- NULL

# Route 1: UNDETECTABLE (highest priority)
if (params$p1_t_mean_range < 0.08 && params$p1_t_mean_range < 0.03) {
  category <- "UNDETECTABLE"
  strategy <- "UNDETECTABLE"
}

# Route 2: PLOIDY_AMBIGUOUS
if (is.null(category) && metrics$p1_multiple_ploidy && 
    metrics$p1_ploidy_loglik_gap < 50) {
  category <- "PLOIDY_AMBIGUOUS"
}

# Route 3: NOISY_WEAK
if (is.null(category) && params$p1_t_mean_range < 0.15 && 
    params$p1_MAD > mad_p90) {
  category <- "NOISY_WEAK"
}

# Route 4: WEAK_SIGNAL
if (is.null(category) && params$p1_t_mean_range < 0.15 && 
    params$p1_MAD <= mad_p90) {
  category <- "WEAK_SIGNAL"
}

# Route 5: NOISY_SIGNAL
if (is.null(category) && params$p1_t_mean_range >= 0.15 && 
    params$p1_MAD > mad_p90) {
  category <- "NOISY_SIGNAL"
}

# Route 6: CLEAN_SIGNAL (default)
if (is.null(category)) {
  category <- "CLEAN_SIGNAL"
}

# ACCEPT overrides all categories
if (p1_accepted) {
  category <- "ACCEPT"
  strategy <- "ACCEPT"
}

# --- Pass 2 parameter construction (skip for ACCEPT/UNDETECTABLE) -----------

pass2_txnE <- NA
pass2_txnStrength <- NA
pass2_maxCN <- NA
pass2_normal <- NA
pass2_ploidy <- NA
pass2_scStates <- NA
pass2_estimateScPrevalence <- NA

if (strategy != "ACCEPT" && strategy != "UNDETECTABLE") {
  # Construct parameters based on category
  param_map <- list(
    UNDETECTABLE = list(
      txnE = 0.9999,
      txnStrength = 10000,
      maxCN = 3,
      scStates = "c()",
      estimateScPrevalence = FALSE
    ),
    NOISY_WEAK = list(
      txnE = 0.9999,
      txnStrength = 10000,
      maxCN = 3,
      scStates = "c()",
      estimateScPrevalence = FALSE
    ),
    WEAK_SIGNAL = list(
      txnE = 0.999,
      txnStrength = 3000,
      maxCN = 4,
      scStates = "c()",
      estimateScPrevalence = FALSE
    ),
    CLEAN_SIGNAL = list(
      txnE = 0.999,
      txnStrength = 3000,
      maxCN = 5,
      scStates = "c()",
      estimateScPrevalence = FALSE
    ),
    NOISY_SIGNAL = list(
      txnE = 0.999,  # Critical: NOT 0.9999
      txnStrength = 5000,
      maxCN = 5,
      scStates = "c()",
      estimateScPrevalence = FALSE
    ),
    PLOIDY_AMBIGUOUS = list(
      txnE = 0.99,
      txnStrength = 1000,
      maxCN = 7,
      scStates = "c(1,3)",
      estimateScPrevalence = TRUE
    )
  )
  
  params_set <- param_map[[category]]
  pass2_txnE <- params_set$txnE
  pass2_txnStrength <- params_set$txnStrength
  pass2_maxCN <- params_set$maxCN
  pass2_scStates <- params_set$scStates
  pass2_estimateScPrevalence <- params_set$estimateScPrevalence
  
  # Construct normal grid
  pass2_normal <- construct_normal_grid(
    params$pass1_tf, params$pass1_ploidy, category, metrics$winner_n_est
  )
  
  # Construct ploidy specification
  if (category == "PLOIDY_AMBIGUOUS") {
    # Lock to dominant (highest loglik) ploidy
    pass2_ploidy <- paste0("c(", round(metrics$winner_phi, 2), ")")
  } else {
    pass2_ploidy <- paste0("c(", round(params$pass1_ploidy, 2), ")")
  }
}

# --- Write strategy.tsv ------------------------------------------------------

strategy_data <- data.frame(
  key = c(
    "strategy",
    "p1_accepted",
    "p1_MAD",
    "p1_t_mean_range",
    "p1_loglik_gap",
    "p1_convergence",
    "p1_at_edge",
    "p1_multiple_ploidy",
    "p1_ploidy_loglik_gap",
    "pass1_tf",
    "pass1_ploidy",
    "pass2_txnE",
    "pass2_txnStrength",
    "pass2_maxCN",
    "pass2_normal",
    "pass2_ploidy",
    "pass2_scStates",
    "pass2_estimateScPrevalence",
    "mad_threshold_source",
    "mad_p90_used"
  ),
  value = c(
    strategy,
    as.character(p1_accepted),
    round(params$p1_MAD, 5),
    round(params$p1_t_mean_range, 4),
    round(metrics$p1_loglik_gap, 2),
    round(metrics$p1_convergence, 4),
    as.character(metrics$p1_at_edge),
    as.character(metrics$p1_multiple_ploidy),
    if (is.finite(metrics$p1_ploidy_loglik_gap)) {
      round(metrics$p1_ploidy_loglik_gap, 2)
    } else {
      "Inf"
    },
    round(params$pass1_tf, 5),
    round(params$pass1_ploidy, 4),
    as.character(pass2_txnE),
    as.character(pass2_txnStrength),
    as.character(pass2_maxCN),
    as.character(pass2_normal),
    as.character(pass2_ploidy),
    as.character(pass2_scStates),
    as.character(pass2_estimateScPrevalence),
    mad_threshold_source,
    round(mad_p90, 5)
  ),
  stringsAsFactors = FALSE
)

# Ensure output directory exists
dir.create(dirname(opt$strategy), recursive = TRUE, showWarnings = FALSE)

# Write TSV
write.table(
  strategy_data,
  opt$strategy,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# --- Console summary ---------------------------------------------------------

`%nin%` <- function(x, y) !(x %in% y)

cat(sprintf("Sample: %s\n", opt$sample))
cat(sprintf("Strategy: %s (category: %s, p1_accepted: %s)\n", 
            strategy, category, p1_accepted))
cat(sprintf("Pass 1 metrics: TF=%.4f, ploidy=%.3f, MAD=%.5f, t_range=%.4f\n",
            params$pass1_tf, params$pass1_ploidy, params$p1_MAD, 
            params$p1_t_mean_range))
cat(sprintf("Loglik gap: %.2f, convergence: %.4f, at_edge: %s\n",
            metrics$p1_loglik_gap, metrics$p1_convergence, metrics$p1_at_edge))
cat(sprintf("MAD threshold: %.5f (source: %s)\n", mad_p90, mad_threshold_source))

if (strategy %nin% c("ACCEPT", "UNDETECTABLE")) {
  cat(sprintf("Pass 2 params: txnE=%s, txnStrength=%s, maxCN=%s\n",
              pass2_txnE, pass2_txnStrength, pass2_maxCN))
  cat(sprintf("  Normal grid: %s\n", pass2_normal))
  cat(sprintf("  Ploidy: %s\n", pass2_ploidy))
}

