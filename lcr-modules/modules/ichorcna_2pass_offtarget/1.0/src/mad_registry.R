#!/usr/bin/env Rscript

# ==============================================================================
# mad_registry.R — MAD Threshold Management and Registry
# ==============================================================================
# Manages the MAD audit trail across batches. Computes P90 and Tukey fence.
# Called as a sourced utility by aggregate_params.R after the params TSV.
# Can also be invoked directly for testing/standalone mode.

# --- Helper: Compute Tukey fence -------------------------------------------------

compute_tukey_fence <- function(values) {
  q1 <- quantile(values, 0.25, names = FALSE)
  q3 <- quantile(values, 0.75, names = FALSE)
  iqr <- q3 - q1
  return(q3 + 1.5 * iqr)
}

# --- Main function: update_mad_registry() ----------------------------------------

update_mad_registry <- function(
  confidence_tsv,
  outdir,
  mad_registry_path = NULL,
  update_registry = FALSE,
  batch_id = NULL,
  run_types = c("off_target", "all_target")
) {
  
  if (!file.exists(confidence_tsv)) {
    stop(sprintf("confidence_tsv not found: %s", confidence_tsv))
  }
  
  # Read current batch confidence TSV (contains p1_MAD column)
  confidence_df <- read.csv(confidence_tsv, stringsAsFactors = FALSE, sep = "\t")
  
  if (!"p1_MAD" %in% names(confidence_df)) {
    stop("confidence_tsv must have a 'p1_MAD' column")
  }
  
  # Infer run_type from filename or strategy (simplified: assume off_target for now)
  if (!"run_type" %in% names(confidence_df)) {
    confidence_df$run_type <- "off_target"  # Default
  }
  
  # Mode selection
  has_registry <- !is.null(mad_registry_path) && file.exists(mad_registry_path)
  
  if (!has_registry) {
    # Mode A: Compute from current batch only
    mode_label <- "A"
    registry_path_used <- NULL
    current_rows <- confidence_df[, c("sample_id", "run_type", "p1_MAD")]
    current_rows$batch_id <- batch_id %||% "current_batch"
    current_rows$date_added <- Sys.Date()
    current_rows$threshold_source <- "current_batch"
    
    # Write registry
    registry_output <- rbind(data.frame(), current_rows)  # Initialize as dataframe
    if (nrow(registry_output) > 0) {
      registry_output$sample_id <- as.character(registry_output$sample_id)
      registry_output$run_type <- as.character(registry_output$run_type)
      registry_output$p1_MAD <- as.numeric(registry_output$p1_MAD)
    }
    
  } else if (!update_registry) {
    # Mode B: Use supplied thresholds, don't recompute
    mode_label <- "B"
    registry_path_used <- mad_registry_path
    
    # Load existing registry
    registry_df <- read.csv(mad_registry_path, stringsAsFactors = FALSE)
    
    # Append current batch rows
    current_rows <- confidence_df[, c("sample_id", "run_type", "p1_MAD")]
    current_rows$batch_id <- batch_id %||% "current_batch"
    current_rows$date_added <- Sys.Date()
    current_rows$threshold_source <- "user_supplied"
    
    # Merge
    registry_output <- rbind(registry_df, current_rows)
    
  } else {
    # Mode C: Merge and recompute
    mode_label <- "C"
    registry_path_used <- mad_registry_path
    
    # Load existing registry
    registry_df <- read.csv(mad_registry_path, stringsAsFactors = FALSE)
    
    # Log pre-merge counts
    pre_counts <- as.data.frame(table(registry_df$run_type))
    names(pre_counts) <- c("run_type", "n_pre")
    
    # Append current batch
    current_rows <- confidence_df[, c("sample_id", "run_type", "p1_MAD")]
    current_rows$batch_id <- batch_id %||% "current_batch"
    current_rows$date_added <- Sys.Date()
    current_rows$threshold_source <- "merged_batch"
    
    registry_output <- rbind(registry_df, current_rows)
    
    # Log post-merge counts
    post_counts <- as.data.frame(table(registry_output$run_type))
    names(post_counts) <- c("run_type", "n_post")
  }
  
  # Compute thresholds per run_type
  thresholds_list <- list()
  
  for (rtype in run_types) {
    run_data <- registry_output[registry_output$run_type == rtype, ]
    
    if (nrow(run_data) == 0) {
      warning(sprintf("No data for run_type=%s", rtype))
      next
    }
    
    mad_values <- run_data$p1_MAD
    mad_values <- mad_values[!is.na(mad_values)]
    
    if (length(mad_values) == 0) {
      warning(sprintf("No valid MAD values for run_type=%s", rtype))
      next
    }
    
    thresholds_list[[rtype]] <- data.frame(
      run_type = rtype,
      n_samples = nrow(run_data),
      mad_median = median(mad_values),
      mad_p75 = quantile(mad_values, 0.75, names = FALSE),
      mad_p90 = quantile(mad_values, 0.90, names = FALSE),
      mad_tukey = compute_tukey_fence(mad_values),
      threshold_source = paste0("Mode_", mode_label),
      registry_path = mad_registry_path %||% "N/A",
      computed_at = Sys.time(),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(thresholds_list) == 0) {
    warning("No thresholds computed")
    return(NULL)
  }
  
  thresholds_df <- do.call(rbind, thresholds_list)
  rownames(thresholds_df) <- NULL
  
  # Write files
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  registry_file <- file.path(outdir, "ichorCNA_mad_registry.csv")
  write.csv(registry_output, registry_file, row.names = FALSE)
  
  thresholds_file <- file.path(outdir, "ichorCNA_mad_thresholds.csv")
  write.csv(thresholds_df, thresholds_file, row.names = FALSE)
  
  # Console logging
  cat(sprintf("Mode %s: MAD thresholds %s\n",
              mode_label,
              if (mode_label == "A") "computed from current batch"
              else if (mode_label == "B") "using pre-computed registry (not recomputed)"
              else "merged and recomputed"))
  cat(sprintf("Registry: %s (%d rows)\n", registry_file, nrow(registry_output)))
  cat(sprintf("Thresholds: %s\n", thresholds_file))
  
  return(list(
    registry = registry_output,
    thresholds = thresholds_df,
    registry_file = registry_file,
    thresholds_file = thresholds_file
  ))
}

# --- Standalone mode (if invoked directly) ----------------------------------------

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    stop("Usage: Rscript mad_registry.R <confidence.tsv> <outdir> [--registry <path>] [--update]")
  }
  
  confidence_file <- args[1]
  output_dir <- args[2]
  
  # Parse optional flags
  reg_path <- NULL
  do_update <- FALSE
  
  if ("--registry" %in% args) {
    idx <- which(args == "--registry")
    if (idx < length(args)) {
      reg_path <- args[idx + 1]
    }
  }
  
  if ("--update" %in% args) {
    do_update <- TRUE
  }
  
  result <- update_mad_registry(
    confidence_file,
    output_dir,
    mad_registry_path = reg_path,
    update_registry = do_update
  )
}

# Helper: %||% operator for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x
