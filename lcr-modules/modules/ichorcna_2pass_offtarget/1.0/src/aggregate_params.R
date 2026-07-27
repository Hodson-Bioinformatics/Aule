#!/usr/bin/env Rscript

# ==============================================================================
# aggregate_params.R — Merge All Per-Sample Solutions Tables
# ==============================================================================
# Reads Pass 1 and Pass 2 params.txt files, merges their solutions tables
# into a single TSV for cohort-level analysis.
#
# Called by rule: _ichorcna_2pass_offtarget_aggregate_params
# CLI: Rscript aggregate_params.R <output.tsv> <p1_params...> <p2_params...>

# --- Helper: Extract sample ID from file path ---------------------------------

extract_sample_id <- function(file_path) {
  # Extract basename, remove .params.txt suffix
  basename <- basename(file_path)
  sample_id <- sub("\\.params\\.txt$", "", basename)
  sample_id <- sub("\\.txt$", "", sample_id)
  return(sample_id)
}

# --- Helper: Parse solutions table from a params file -------------------------

parse_solutions_from_params <- function(params_file, sample_id, pass_num) {
  if (!file.exists(params_file)) {
    stop(sprintf("params file not found: %s", params_file))
  }
  
  lines <- readLines(params_file)
  
  # Find the "init\tn_est\tphi_est\t..." header line
  init_idx <- grep("^init\\b", lines)[1]
  if (is.na(init_idx)) {
    warning(sprintf("Could not find solutions table in %s", params_file))
    return(NULL)
  }
  
  # Extract solution rows
  init_lines <- lines[(init_idx + 1):length(lines)]
  init_lines <- init_lines[nchar(trimws(init_lines)) > 0]
  
  # Parse each row
  solutions <- list()
  for (line in init_lines) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) < 7) next
    
    # Parse columns: init, n_est, phi_est, BIC, Frac_genome_subclonal, Frac_CNA_subclonal, loglik
    loglik_val <- tryCatch(as.numeric(parts[7]), error = function(e) NA)
    if (is.na(loglik_val)) next
    
    solutions[[length(solutions) + 1]] <- data.frame(
      init = parts[1],
      n_est = as.numeric(parts[2]),
      phi_est = as.numeric(parts[3]),
      BIC = if (is.na(as.numeric(parts[4]))) NA else as.numeric(parts[4]),
      Frac_genome_subclonal = if (is.na(as.numeric(parts[5]))) NA else as.numeric(parts[5]),
      Frac_CNA_subclonal = if (is.na(as.numeric(parts[6]))) NA else as.numeric(parts[6]),
      loglik = loglik_val,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(solutions) == 0) {
    warning(sprintf("No solutions parsed from %s", params_file))
    return(NULL)
  }
  
  solutions_df <- do.call(rbind, solutions)
  rownames(solutions_df) <- NULL
  
  # Add metadata columns
  solutions_df$sample_id <- sample_id
  solutions_df$pass <- pass_num
  
  # Determine if this is the winner (highest loglik)
  max_loglik <- max(solutions_df$loglik)
  solutions_df$is_winner <- (solutions_df$loglik == max_loglik)
  
  return(solutions_df)
}

# --- Main execution ----------------------------------------------------------

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript aggregate_params.R <output.tsv> <p1_params...> <p2_params...>")
}

output_path <- args[1]
all_params_files <- args[-1]  # All remaining arguments

# Determine the split point: equal number of P1 and P2 files
n_params <- length(all_params_files)
if (n_params %% 2 != 0) {
  stop("Error: expected equal number of Pass 1 and Pass 2 params files")
}

n_samples <- n_params / 2
p1_files <- all_params_files[1:n_samples]
p2_files <- all_params_files[(n_samples + 1):(2 * n_samples)]

# Parse all solutions tables
all_solutions <- list()

for (i in seq_along(p1_files)) {
  sample_id <- extract_sample_id(p1_files[i])
  
  # Parse Pass 1
  p1_solutions <- parse_solutions_from_params(p1_files[i], sample_id, 1)
  if (!is.null(p1_solutions)) {
    all_solutions[[length(all_solutions) + 1]] <- p1_solutions
  }
  
  # Parse Pass 2
  p2_solutions <- parse_solutions_from_params(p2_files[i], sample_id, 2)
  if (!is.null(p2_solutions)) {
    all_solutions[[length(all_solutions) + 1]] <- p2_solutions
  }
}

# Merge all solutions into a single dataframe
if (length(all_solutions) == 0) {
  stop("No solutions were successfully parsed")
}

merged_df <- do.call(rbind, all_solutions)
rownames(merged_df) <- NULL

# Reorder columns for clarity
col_order <- c(
  "sample_id", "pass", "init", "n_est", "phi_est", 
  "BIC", "Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik", "is_winner"
)

# Only include columns that exist
col_order <- col_order[col_order %in% names(merged_df)]
merged_df <- merged_df[, col_order]

# Ensure output directory exists
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Write output
write.table(
  merged_df,
  output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Console summary
cat(sprintf("Merged %d samples × 2 passes = %d rows\n", n_samples, nrow(merged_df)))
cat(sprintf("Output: %s\n", output_path))
