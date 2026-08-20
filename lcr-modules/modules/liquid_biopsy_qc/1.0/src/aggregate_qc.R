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

safe_spearman <- function(x, y) {
  x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(c(n = sum(keep), rho = NA_real_))
  c(n = sum(keep), rho = unname(cor(x[keep], y[keep], method = "spearman")))
}

burden_correlation_table <- function(table, threshold = 0.30) {
  targets <- c("burden_amp", "tf_estimate")
  missing <- setdiff(targets, names(table))
  if (length(missing)) stop("Missing correlation target(s): ", paste(missing, collapse = ", "), call. = FALSE)
  numeric_names <- names(table)[vapply(table, function(x) {
    values <- suppressWarnings(as.numeric(as.character(x)))
    sum(is.finite(values)) > 0L
  }, logical(1))]
  numeric_names <- setdiff(numeric_names, targets)
  rows <- lapply(numeric_names, function(metric) {
    against_burden <- safe_spearman(table[[metric]], table$burden_amp)
    against_tf <- safe_spearman(table[[metric]], table$tf_estimate)
    confounded <- any(abs(c(against_burden["rho"], against_tf["rho"])) > threshold,
                      na.rm = TRUE) || identical(metric, "restart_convergence")
    data.frame(
      metric = metric,
      n_burden_amp = as.integer(against_burden["n"]),
      rho_burden_amp = against_burden["rho"],
      n_tf_estimate = as.integer(against_tf["n"]),
      rho_tf_estimate = against_tf["rho"],
      verdict = if (confounded) {
        "burden-confounded - not eligible as a filter"
      } else {
        "not burden-confounded at configured threshold"
      },
      filter_eligible = !confounded,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

cohort_summary <- function(table) {
  numeric_names <- names(table)[vapply(table, is.numeric, logical(1))]
  rows <- lapply(numeric_names, function(metric) {
    values <- table[[metric]][is.finite(table[[metric]])]
    q <- if (length(values)) quantile(values, c(0, .25, .5, .75, .9, 1), names = FALSE) else rep(NA_real_, 6)
    data.frame(summary_type = "metric", metric = metric, n = length(values),
               q0 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q90 = q[5], q100 = q[6],
               value = NA_character_, stringsAsFactors = FALSE)
  })
  flags <- unlist(strsplit(as.character(table$qc_flag), ";", fixed = TRUE))
  counts <- sort(base::table(flags), decreasing = TRUE)
  flag_rows <- data.frame(summary_type = "qc_flag", metric = names(counts), n = as.integer(counts),
                          q0 = NA_real_, q25 = NA_real_, q50 = NA_real_, q75 = NA_real_,
                          q90 = NA_real_, q100 = NA_real_, value = NA_character_,
                          stringsAsFactors = FALSE)
  bin_values <- unique(table$n_bins_used[is.finite(table$n_bins_used)])
  bin_row <- data.frame(summary_type = "assertion", metric = "n_bins_used_constant",
                        n = sum(is.finite(table$n_bins_used)), q0 = NA_real_, q25 = NA_real_,
                        q50 = NA_real_, q75 = NA_real_, q90 = NA_real_, q100 = NA_real_,
                        value = if (length(bin_values) <= 1L) "pass" else "FAIL",
                        stringsAsFactors = FALSE)
  do.call(rbind, c(rows, list(flag_rows, bin_row)))
}

with_device <- function(pdf_path, png_path, draw) {
  width <- 89 / 25.4; height <- 70 / 25.4
  pdf(pdf_path, width = width, height = height, family = "sans", pointsize = 7.5,
      useDingbats = FALSE)
  par(mar = c(3.2, 3.3, 1.0, 0.5), mgp = c(1.8, .5, 0), tcl = -.2,
      las = 1, cex.axis = .85, cex.lab = .9)
  draw(); dev.off()
  png(png_path, width = round(width * 600), height = round(height * 600),
      res = 600, type = "cairo", pointsize = 7.5)
  par(mar = c(3.2, 3.3, 1.0, 0.5), mgp = c(1.8, .5, 0), tcl = -.2,
      las = 1, cex.axis = .85, cex.lab = .9)
  draw(); dev.off()
}

blank_plot <- function(label) {
  plot.new(); text(.5, .5, label, cex = .9)
}

read_insert_histograms <- function(raw_dir) {
  paths <- list.files(raw_dir, pattern = "[.]insert_size_metrics[.]txt$",
                      recursive = TRUE, full.names = TRUE)
  output <- list()
  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    header <- grep("^insert_size\\t", lines, ignore.case = TRUE)[1L]
    if (is.na(header)) next
    rows <- lines[(header + 1L):length(lines)]
    rows <- rows[grepl("^[0-9]+\\t", rows)]
    if (!length(rows)) next
    tab <- read.table(text = paste(c(lines[header], rows), collapse = "\n"),
                      sep = "\t", header = TRUE, check.names = FALSE)
    output[[sub("[.]insert_size_metrics[.]txt$", "", basename(path))]] <-
      data.frame(size = as.numeric(tab[[1L]]), count = as.numeric(tab[[2L]]))
  }
  output
}

draw_figures <- function(table, raw_dir, pon_reference, figure_dir) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"; grey <- "#77777755"
  path <- function(stem, ext) file.path(figure_dir, paste0(stem, ".", ext))

  histograms <- read_insert_histograms(raw_dir)
  with_device(path("insert_size_distributions", "pdf"), path("insert_size_distributions", "png"), function() {
    if (!length(histograms)) return(blank_plot("Insert-size metrics unavailable"))
    xmax <- max(unlist(lapply(histograms, function(x) x$size[x$size <= 500])), na.rm = TRUE)
    plot(NA, xlim = c(0, xmax), ylim = c(0, 1), xlab = "Insert size (bp)", ylab = "Relative density")
    for (hist in histograms) {
      keep <- hist$size <= xmax & is.finite(hist$count)
      lines(hist$size[keep], hist$count[keep] / max(hist$count[keep]), col = grey, lwd = .35)
    }
    abline(v = 166, col = orange, lty = 2, lwd = .7)
  })

  reference <- tryCatch(read_tsv(pon_reference), error = function(e) NULL)
  ref_noise <- if (!is.null(reference) && "sigma_noise" %in% names(reference)) {
    reference$sigma_noise[is.finite(reference$sigma_noise)]
  } else numeric()
  with_device(path("sigma_noise_vs_pon", "pdf"), path("sigma_noise_vs_pon", "png"), function() {
    y <- table$sigma_noise
    if (!any(is.finite(y))) return(blank_plot("Noise estimates unavailable"))
    plot(seq_along(y), y, pch = 16, cex = .55, col = blue,
         xlab = "Sample (ordered)", ylab = expression(sigma[noise]))
    if (length(ref_noise)) rect(par("usr")[1], min(ref_noise), par("usr")[2], max(ref_noise),
                                col = "#009E7330", border = NA)
    points(seq_along(y), y, pch = 16, cex = .55, col = blue)
  })

  with_device(path("resid_gc_cor_vs_burden_amp", "pdf"), path("resid_gc_cor_vs_burden_amp", "png"), function() {
    keep <- is.finite(table$burden_amp) & is.finite(table$resid_gc_cor)
    if (sum(keep) < 2L) return(blank_plot("Residual-GC metrics unavailable"))
    plot(table$burden_amp[keep], table$resid_gc_cor[keep], pch = 16, cex = .55,
         col = blue, xlab = "Burden amplitude", ylab = "Residual GC Spearman rho")
    abline(h = 0, col = grey, lwd = .5)
  })

  flags <- unlist(strsplit(as.character(table$qc_flag), ";", fixed = TRUE))
  flags <- flags[flags != "pass"]
  counts <- sort(base::table(flags), decreasing = TRUE)
  with_device(path("flag_frequency", "pdf"), path("flag_frequency", "png"), function() {
    if (!length(counts)) return(blank_plot("No advisory QC flags"))
    par(mar = c(3.2, 7, 1, .5))
    barplot(rev(counts), horiz = TRUE, las = 1, col = green, border = NA,
            xlab = "Samples", cex.names = .75)
  })
}

aggregate_qc <- function(paths, rho_threshold, raw_dir, pon_reference, figure_dir) {
  tables <- lapply(paths, read_tsv)
  all_names <- Reduce(union, lapply(tables, names))
  fill <- function(table) {
    for (name in setdiff(all_names, names(table))) table[[name]] <- NA
    table[, all_names, drop = FALSE]
  }
  cohort <- do.call(rbind, lapply(tables, fill))
  rownames(cohort) <- NULL
  numeric_candidates <- setdiff(names(cohort), c(
    "sample_id", "assay_mode", "coverage_scope", "qc_flag", "missing_reason", "pon_id", "pon_md5",
    "gc_wig_md5", "map_wig_md5", "ichorcna_commit_sha", "genome_build",
    "genomeStyle", "pon_panel", "pon_library_prep", "pon_cohort"
  ))
  for (name in numeric_candidates) {
    converted <- suppressWarnings(as.numeric(as.character(cohort[[name]])))
    if (sum(is.finite(converted)) || all(is.na(cohort[[name]]))) cohort[[name]] <- converted
  }
  correlations <- burden_correlation_table(cohort, rho_threshold)
  summary <- cohort_summary(cohort)
  provenance_names <- intersect(c("sample_id", "assay_mode", "coverage_scope", "pon_id", "pon_md5", "pon_n_normals",
                                  "gc_wig_md5", "map_wig_md5", "ichorcna_commit_sha",
                                  "binsize_kb", "genome_build", "genomeStyle", "pon_panel",
                                  "pon_library_prep", "pon_cohort"), names(cohort))
  draw_figures(cohort, raw_dir, pon_reference, figure_dir)
  list(cohort = cohort, summary = summary, correlations = correlations,
       provenance = cohort[, provenance_names, drop = FALSE])
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_args(args)
  paths <- strsplit(options$inputs, ",", fixed = TRUE)[[1L]]
  result <- aggregate_qc(paths, as.numeric(options$`rho-threshold`), options$`raw-dir`,
                         options$`pon-reference`, options$`figure-dir`)
  for (path in c(options$`per-sample`, options$summary, options$correlation, options$provenance))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(result$cohort, options$`per-sample`, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(result$summary, options$summary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(result$correlations, options$correlation, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(result$provenance, options$provenance, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

if (!interactive() && sys.nframe() == 0L && length(commandArgs(trailingOnly = TRUE))) main()
