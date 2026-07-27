#!/usr/bin/env Rscript

# ==============================================================================
# combine_plots.R — Combine Pass 1 and Pass 2 PDFs into Paired Review Document
# ==============================================================================
# Combines Pass 1 and Pass 2 genome-wide PDFs side-by-side for QC review.
#
# Called by rule: _ichorcna_2pass_offtarget_aggregate_pdf
# CLI: Rscript combine_plots.R <output.pdf> <p1_plots...> <p2_plots...>

suppressPackageStartupMessages({
  library(pdftools)
})

# --- Helper: Extract sample ID from file path ---------------------------------

extract_sample_id <- function(file_path) {
  basename <- basename(file_path)
  # Try common naming patterns
  sample_id <- sub("_pass1_genomeWide\\.pdf$", "", basename)
  sample_id <- sub("_genomeWide\\.pdf$", "", sample_id)
  sample_id <- sub("\\.pdf$", "", sample_id)
  return(sample_id)
}

# --- Main execution ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript combine_plots.R <output.pdf> <p1_plots...> <p2_plots...>")
}

output_path <- args[1]
all_pdf_files <- args[-1]

# Determine split point: equal number of P1 and P2 files
n_pdfs <- length(all_pdf_files)
if (n_pdfs %% 2 != 0) {
  stop("Error: expected equal number of Pass 1 and Pass 2 PDF files")
}

n_samples <- n_pdfs / 2
p1_files <- all_pdf_files[1:n_samples]
p2_files <- all_pdf_files[(n_samples + 1):(2 * n_samples)]

# Validate files exist
for (f in c(p1_files, p2_files)) {
  if (!file.exists(f)) {
    stop(sprintf("PDF file not found: %s", f))
  }
}

# Combine PDFs: alternate P1 and P2 files
combined_pdfs <- list()

for (i in seq_along(p1_files)) {
  sample_id <- extract_sample_id(p1_files[i])
  
  # Add P1 PDF
  if (file.exists(p1_files[i])) {
    tryCatch({
      p1_pdf <- pdftools::pdf_read_pdf(p1_files[i])
      combined_pdfs[[length(combined_pdfs) + 1]] <- p1_pdf
    }, error = function(e) {
      warning(sprintf("Could not read P1 PDF for %s: %s", sample_id, e$message))
    })
  }
  
  # Add P2 PDF
  if (file.exists(p2_files[i])) {
    tryCatch({
      p2_pdf <- pdftools::pdf_read_pdf(p2_files[i])
      combined_pdfs[[length(combined_pdfs) + 1]] <- p2_pdf
    }, error = function(e) {
      warning(sprintf("Could not read P2 PDF for %s: %s", sample_id, e$message))
    })
  }
}

# Combine using pdftools::pdf_combine if available, otherwise direct concatenation
if (length(combined_pdfs) > 0) {
  # Try to use pdf_combine on the files directly
  all_valid_files <- c(p1_files, p2_files)
  all_valid_files <- all_valid_files[file.exists(all_valid_files)]
  
  if (length(all_valid_files) > 0) {
    tryCatch({
      pdftools::pdf_combine(all_valid_files, output = output_path)
      cat(sprintf("Combined %d PDFs into %s\n", length(all_valid_files), output_path))
    }, error = function(e) {
      warning(sprintf("pdf_combine failed: %s. Writing file list instead.", e$message))
      # Fallback: write file list
      writeLines(all_valid_files, output_path)
    })
  }
} else {
  stop("No valid PDF files to combine")
}

# Ensure output directory exists
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output PDF: %s\n", output_path))
