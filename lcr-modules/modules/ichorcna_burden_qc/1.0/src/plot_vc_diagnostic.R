parse_plot_args <- function(args) {
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
      if (index == length(args) || startsWith(args[index + 1L], "--")) {
        stop("Missing value for --", key, call. = FALSE)
      }
      index <- index + 1L
      value <- args[index]
    }
    parsed[[key]] <- value
    index <- index + 1L
  }
  parsed
}

read_vc_diagnostic <- function(table_path, sample_id) {
  if (!file.exists(table_path)) {
    stop("Burden table does not exist: ", table_path, call. = FALSE)
  }
  burden <- read.table(
    table_path, sep = "\t", header = TRUE, check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE,
    colClasses = "character"
  )
  required <- c(
    "sample_id", "window_widths", "burden_var_signed", "burden_amp",
    "sigma_noise", "vc_fit_r2"
  )
  missing_required <- setdiff(required, names(burden))
  if (length(missing_required) > 0L) {
    stop(
      "Burden table is missing required column(s): ",
      paste(missing_required, collapse = ", "), call. = FALSE
    )
  }
  matches <- which(as.character(burden$sample_id) == sample_id)
  if (length(matches) == 0L) {
    stop("sample_id not found in burden table: ", sample_id, call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop("sample_id matches multiple rows in burden table: ", sample_id,
         call. = FALSE)
  }
  row <- burden[matches, , drop = FALSE]
  width_text <- trimws(strsplit(
    as.character(row$window_widths), ",", fixed = TRUE
  )[[1L]])
  widths <- suppressWarnings(as.numeric(width_text))
  if (length(widths) == 0L || any(!is.finite(widths)) ||
      any(widths <= 0) || anyDuplicated(widths)) {
    stop("Invalid window_widths for sample ", sample_id, ": ",
         row$window_widths, call. = FALSE)
  }
  width_suffixes <- as.character(widths)
  variance_columns <- paste0("vc_V_w", width_suffixes)
  count_columns <- paste0("vc_nwin_w", width_suffixes)
  missing_variances <- setdiff(variance_columns, names(row))
  if (length(missing_variances) > 0L) {
    stop(
      "Variance-component columns are absent (old burden TSV). Missing: ",
      paste(missing_variances, collapse = ", "), call. = FALSE
    )
  }
  missing_counts <- setdiff(count_columns, names(row))
  if (length(missing_counts) > 0L) {
    stop(
      "Variance-component window-count columns are missing: ",
      paste(missing_counts, collapse = ", "), call. = FALSE
    )
  }
  numeric_value <- function(name) {
    suppressWarnings(as.numeric(row[[name]]))
  }
  variances <- vapply(variance_columns, numeric_value, numeric(1))
  counts <- vapply(count_columns, numeric_value, numeric(1))
  intercept <- numeric_value("burden_var_signed")
  sigma_noise <- numeric_value("sigma_noise")
  burden_amp <- numeric_value("burden_amp")
  fit_r2 <- numeric_value("vc_fit_r2")
  if (!any(is.finite(variances))) {
    stop("No finite vc_V_w* values for sample: ", sample_id, call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("vc_nwin_w* values must be finite and non-negative for sample: ",
         sample_id, call. = FALSE)
  }
  if (!is.finite(intercept) || !is.finite(sigma_noise) ||
      !is.finite(burden_amp) || !is.finite(fit_r2)) {
    stop("Stored burden fit values are non-finite for sample: ", sample_id,
         call. = FALSE)
  }
  list(
    sample_id = sample_id,
    widths = widths,
    variances = variances,
    n_win = counts,
    burden_var_signed = intercept,
    burden_amp = burden_amp,
    sigma_noise = sigma_noise,
    vc_fit_r2 = fit_r2
  )
}

open_plot_device <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (!extension %in% c("png", "pdf")) {
    stop("--out must end in .png or .pdf", call. = FALSE)
  }
  output_directory <- dirname(path)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  if (extension == "png") {
    grDevices::png(path, width = 1600, height = 1200, res = 200,
                   bg = "white")
  } else {
    grDevices::pdf(path, width = 8, height = 6, onefile = FALSE)
  }
}

plot_vc_diagnostic <- function(table_path, sample_id, out_path) {
  diagnostic <- read_vc_diagnostic(table_path, sample_id)
  open_plot_device(out_path)
  on.exit(grDevices::dev.off(), add = TRUE)

  inverse_widths <- 1 / diagnostic$widths
  finite_points <- is.finite(diagnostic$variances)
  max_variance <- max(diagnostic$variances[finite_points])
  y_upper <- if (max_variance > 0) max_variance * 1.08 else 1
  root_counts <- sqrt(diagnostic$n_win)
  point_cex <- if (max(root_counts) > 0) {
    2.4 * root_counts / max(root_counts)
  } else {
    rep(1, length(root_counts))
  }
  slope <- diagnostic$sigma_noise^2

  graphics::plot(
    NA_real_, NA_real_, xlim = c(0, 1.04), ylim = c(0, y_upper),
    xaxs = "i", yaxs = "i", xlab = "Inverse window width (1/w)",
    ylab = "Variance of window means", axes = FALSE
  )
  graphics::axis(1, at = sort(unique(c(0, inverse_widths))))
  graphics::axis(2, las = 1)
  graphics::box()
  graphics::segments(
    0, diagnostic$burden_var_signed,
    1, diagnostic$burden_var_signed + slope,
    col = "#0072B2", lwd = 2
  )
  graphics::points(
    inverse_widths[finite_points], diagnostic$variances[finite_points],
    pch = 21, bg = "#E69F00", col = "black",
    cex = point_cex[finite_points], lwd = 1
  )
  graphics::points(
    0, diagnostic$burden_var_signed, pch = 23,
    bg = "#CC79A7", col = "black", cex = 1.5, lwd = 1.2
  )
  graphics::text(
    0.025, diagnostic$burden_var_signed,
    labels = sprintf("burden_amp = sqrt(S) = %.4g", diagnostic$burden_amp),
    pos = 4, xpd = NA
  )
  annotation_x <- 0.68
  slope_label <- bquote(
    sigma^2~"(technical noise)" == .(format(slope, digits = 4))
  )
  graphics::text(
    annotation_x,
    diagnostic$burden_var_signed + slope * annotation_x,
    labels = slope_label,
    pos = 3, col = "#0072B2"
  )
  graphics::title(main = diagnostic$sample_id, line = 2)
  graphics::mtext(
    sprintf("vc_fit_r2 = %.4f", diagnostic$vc_fit_r2),
    side = 3, line = 0.35
  )

  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  invisible(diagnostic)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_plot_args(args)
  required <- c("table", "sample", "out")
  missing <- required[!required %in% names(options)]
  if (length(missing) > 0L) {
    stop("Missing required argument(s): ",
         paste(paste0("--", missing), collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(names(options), required)
  if (length(unknown) > 0L) {
    stop("Unknown argument(s): ",
         paste(paste0("--", unknown), collapse = ", "), call. = FALSE)
  }
  plot_vc_diagnostic(options$table, options$sample, options$out)
}

if (!interactive() && sys.nframe() == 0L &&
    length(commandArgs(trailingOnly = TRUE)) > 0L) main()
