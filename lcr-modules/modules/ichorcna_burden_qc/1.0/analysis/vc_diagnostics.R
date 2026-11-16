suppressPackageStartupMessages(library(ggplot2))

VC_DIAGNOSTIC_DATA <- NULL
VC_DIAGNOSTIC_LONG <- NULL

parse_args <- function(args) {
  output <- list()
  index <- 1L
  while (index <= length(args)) {
    key <- sub("^--", "", args[index])
    index <- index + 1L
    output[[key]] <- args[index]
    index <- index + 1L
  }
  output
}

read_tsv <- function(path) {
  read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
             quote = "", comment.char = "", stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
}

backfill_vc_diagnostics <- function(burden, corrected_dir) {
  paths <- list.files(corrected_dir, pattern = "[.]correctedDepth[.]txt$",
                      recursive = TRUE, full.names = TRUE)
  ids <- sub("[.]correctedDepth[.]txt$", "", basename(paths))
  path_by_id <- setNames(paths, ids)
  diagnostic_rows <- lapply(seq_len(nrow(burden)), function(index) {
    sample_id <- burden$sample_id[index]
    prepared <- prepare_bins(path_by_id[[sample_id]])
    widths <- as.numeric(strsplit(
      burden$window_widths[index], ",", fixed = TRUE
    )[[1L]])
    amplitude <- burden_amplitude(
      prepared$bins$chr, prepared$bins$logr, widths,
      burden$winsor_p[index]
    )
    cbind(
      data.frame(
        sample_id = sample_id,
        vc_curv = amplitude$vc_curv,
        S_fine = amplitude$S_fine,
        S_coarse = amplitude$S_coarse,
        S_excess = amplitude$S_excess,
        S_unweighted = amplitude$S_unweighted,
        amp_coarse = amplitude$amp_coarse,
        stringsAsFactors = FALSE
      ),
      vc_diagnostic_columns(amplitude)
    )
  })
  diagnostics <- do.call(rbind, diagnostic_rows)
  replace_names <- setdiff(names(diagnostics), "sample_id")
  burden[replace_names] <- NULL
  merge(burden, diagnostics, by = "sample_id", sort = FALSE)
}

make_long_data <- function(data) {
  value_names <- grep("^vc_V_w", names(data), value = TRUE)
  rows <- lapply(seq_len(nrow(data)), function(index) {
    widths <- as.numeric(sub("vc_V_w", "", value_names))
    values <- as.numeric(data[index, value_names])
    counts <- as.numeric(data[index, paste0("vc_nwin_w", widths)])
    residuals <- as.numeric(data[index, paste0("vc_resid_w", widths)])
    v1 <- values[widths == 1]
    fitted <- data$S_fine[index] + data$sigma_noise[index]^2 / widths
    data.frame(
      sample_id = data$sample_id[index],
      burden_amp = data$burden_amp[index],
      sigma_noise = data$sigma_noise[index],
      width = widths,
      inv_w = 1 / widths,
      variance = values,
      n_win = counts,
      residual = residuals,
      variance_norm = values / v1,
      fitted = fitted,
      fitted_norm = fitted / v1,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

nature_theme <- function() {
  theme_classic(base_family = "sans", base_size = 7.5) +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.5, "mm"),
      legend.key.height = unit(4, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 7.5),
      plot.title = element_text(size = 8, face = "bold"),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

save_figure <- function(plot, stem, width_mm, height_mm) {
  ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(stem, ".png"), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, bg = "white")
}

plot_cohort_overlay <- function(long) {
  mean_data <- aggregate(
    cbind(variance_norm, fitted_norm) ~ width + inv_w,
    long, mean, na.rm = TRUE
  )
  ggplot(long, aes(inv_w, variance_norm, group = sample_id,
                   colour = sigma_noise)) +
    geom_line(linewidth = 0.2, alpha = 0.16) +
    geom_line(data = mean_data, aes(inv_w, variance_norm, group = 1),
              inherit.aes = FALSE, colour = "black", linewidth = 0.8) +
    geom_line(data = mean_data, aes(inv_w, fitted_norm, group = 1),
              inherit.aes = FALSE, colour = "black", linewidth = 0.55,
              linetype = "22") +
    scale_colour_viridis_c(option = "D", end = 0.92,
                           name = expression(sigma[noise])) +
    scale_x_continuous(breaks = 1 / c(1, 2, 4, 8, 16),
                       labels = c("1", "1/2", "1/4", "1/8", "1/16")) +
    labs(x = "Inverse window width (1/w)", y = "V(w) / V(1)") +
    nature_theme()
}

residual_summary <- function(long) {
  split_data <- split(long, long$width)
  do.call(rbind, lapply(split_data, function(group) {
    values <- group$residual[is.finite(group$residual)]
    standard_error <- sd(values) / sqrt(length(values))
    data.frame(
      width = group$width[1L],
      inv_w = group$inv_w[1L],
      mean_residual = mean(values),
      lo95 = mean(values) - qt(0.975, length(values) - 1L) * standard_error,
      hi95 = mean(values) + qt(0.975, length(values) - 1L) * standard_error,
      n = length(values)
    )
  }))
}

plot_residual_structure <- function(summary) {
  ggplot(summary, aes(inv_w, mean_residual)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35,
               linetype = "22") +
    geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.015,
                  linewidth = 0.35) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.4, shape = 21, fill = "white", stroke = 0.35) +
    scale_x_continuous(breaks = 1 / c(1, 2, 4, 8, 16),
                       labels = c("1", "1/2", "1/4", "1/8", "1/16")) +
    labs(x = "Inverse window width (1/w)", y = "Mean fit residual") +
    nature_theme()
}

choose_exemplars <- function(data) {
  ordered <- data[order(data$burden_amp), ]
  groups <- split(seq_len(nrow(ordered)),
                  cut(seq_len(nrow(ordered)), 3, labels = FALSE))
  burden_labels <- c("Low burden", "Median burden", "High burden")
  selected <- do.call(rbind, lapply(seq_along(groups), function(index) {
    group <- ordered[groups[[index]], ]
    low <- group[which.min(group$sigma_noise), ]
    high <- group[which.max(group$sigma_noise), ]
    low$burden_band <- high$burden_band <- burden_labels[index]
    low$noise_band <- "Low noise"
    high$noise_band <- "High noise"
    rbind(low, high)
  }))
  selected$burden_band <- factor(selected$burden_band, levels = burden_labels)
  selected$noise_band <- factor(selected$noise_band,
                                levels = c("Low noise", "High noise"))
  selected
}

single_sample_plot <- function(sample_id) {
  sample_data <- VC_DIAGNOSTIC_DATA[
    VC_DIAGNOSTIC_DATA$sample_id == sample_id, , drop = FALSE
  ]
  points <- VC_DIAGNOSTIC_LONG[
    VC_DIAGNOSTIC_LONG$sample_id == sample_id, , drop = FALSE
  ]
  fit <- data.frame(
    inv_w = c(0, 1),
    variance = sample_data$S_fine + sample_data$sigma_noise^2 * c(0, 1)
  )
  ggplot(points, aes(inv_w, variance)) +
    geom_line(data = fit, linewidth = 0.45, colour = "#0072B2") +
    geom_point(aes(size = n_win), shape = 21, fill = "#E69F00",
               stroke = 0.35) +
    geom_point(data = data.frame(inv_w = 0, variance = sample_data$S_fine),
               aes(x = inv_w, y = variance), shape = 4, size = 2,
               stroke = 0.55, inherit.aes = FALSE) +
    annotate(
      "text", x = 0.03, y = Inf, hjust = 0, vjust = 1.4, size = 2.3,
      label = sprintf("burden_amp = %.3f", sample_data$burden_amp)
    ) +
    scale_size_continuous(range = c(1.1, 3.2), guide = "none") +
    scale_x_continuous(breaks = c(0, 1 / c(16, 8, 4, 2, 1))) +
    labs(title = sample_id, x = "Inverse window width (1/w)",
         y = "Variance of window means") +
    nature_theme()
}

plot_vc_diagnostic <- function(sample_id, out_path) {
  plot <- single_sample_plot(sample_id)
  extension <- tolower(tools::file_ext(out_path))
  if (extension == "pdf") {
    ggsave(out_path, plot, width = 89, height = 72, units = "mm",
           device = cairo_pdf)
  } else {
    ggsave(out_path, plot, width = 89, height = 72, units = "mm",
           dpi = 600, bg = "white")
  }
  invisible(plot)
}

plot_exemplars <- function(data, long) {
  exemplars <- choose_exemplars(data)
  panel_data <- merge(
    long, exemplars[c("sample_id", "burden_band", "noise_band")],
    by = "sample_id"
  )
  fits <- do.call(rbind, lapply(seq_len(nrow(exemplars)), function(index) {
    sample <- exemplars[index, ]
    data.frame(
      sample_id = sample$sample_id,
      burden_band = sample$burden_band,
      noise_band = sample$noise_band,
      inv_w = c(0, 1),
      variance = sample$S_fine + sample$sigma_noise^2 * c(0, 1)
    )
  }))
  annotations <- transform(
    exemplars, inv_w = 0.03, variance = Inf,
    label = sprintf("%s\namp = %.3f", sample_id, burden_amp)
  )
  plot <- ggplot(panel_data, aes(inv_w, variance)) +
    geom_line(data = fits, linewidth = 0.4, colour = "#0072B2") +
    geom_point(aes(size = n_win), shape = 21, fill = "#E69F00",
               stroke = 0.3) +
    geom_point(data = exemplars,
               aes(x = 0, y = S_fine), inherit.aes = FALSE,
               shape = 4, size = 1.7, stroke = 0.45) +
    geom_text(data = annotations,
              aes(inv_w, variance, label = label), inherit.aes = FALSE,
              hjust = 0, vjust = 1.25, size = 2.15) +
    facet_grid(noise_band ~ burden_band, scales = "free_y") +
    scale_size_continuous(range = c(0.8, 2.6), guide = "none") +
    labs(x = "Inverse window width (1/w)", y = "Variance of window means") +
    nature_theme()
  list(plot = plot, exemplars = exemplars)
}

plot_mechanism <- function(data, clinical) {
  merged <- merge(data, clinical, by = "sample_id")
  threshold <- quantile(merged$PVbest_ctDNA, 0.2, na.rm = TRUE)
  low <- merged[merged$PVbest_ctDNA <= threshold, ]
  ggplot(merged, aes(sigma_noise^2, S_excess, colour = PVbest_ctDNA)) +
    geom_point(size = 1.2, alpha = 0.75) +
    geom_smooth(data = low, method = "lm", formula = y ~ x + 0,
                se = TRUE, colour = "black", fill = "grey75",
                linewidth = 0.5) +
    scale_colour_viridis_c(option = "D", name = "PVbest ctDNA") +
    labs(x = expression(sigma[noise]^2), y = expression(S[fine] - S[coarse])) +
    nature_theme()
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_args(args)
  script_path <- normalizePath(sub("^--file=", "", grep(
    "^--file=", commandArgs(), value = TRUE
  )[1L]))
  source(file.path(dirname(script_path), "../src/burden_index.R"))
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)

  data <- read_tsv(options$burden)
  if (!all(c("vc_V_w1", "S_coarse", "vc_curv") %in% names(data))) {
    data <- backfill_vc_diagnostics(data, options$corrected_dir)
  }
  long <- make_long_data(data)
  VC_DIAGNOSTIC_DATA <<- data
  VC_DIAGNOSTIC_LONG <<- long

  write_tsv(data, file.path(options$output_dir, "vc_diagnostics.tsv"))
  write_tsv(long, file.path(options$output_dir, "vc_window_data.tsv"))

  overlay <- plot_cohort_overlay(long)
  save_figure(overlay, file.path(options$output_dir, "figure1_cohort_overlay"),
              183, 105)

  residuals <- residual_summary(long)
  write_tsv(residuals, file.path(options$output_dir, "vc_residual_summary.tsv"))
  residual_plot <- plot_residual_structure(residuals)
  save_figure(residual_plot,
              file.path(options$output_dir, "figure2_residual_structure"),
              89, 75)

  exemplar_result <- plot_exemplars(data, long)
  write_tsv(exemplar_result$exemplars,
            file.path(options$output_dir, "vc_exemplars.tsv"))
  save_figure(exemplar_result$plot,
              file.path(options$output_dir, "figure3_exemplars"), 183, 125)

  if (!is.null(options$clinical)) {
    mechanism <- plot_mechanism(data, read_tsv(options$clinical))
    save_figure(mechanism,
                file.path(options$output_dir, "figure4_mechanism"), 89, 78)
  }

  sink(file.path(options$output_dir, "session_info.txt"))
  cat("seed: 1\n")
  print(sessionInfo())
  sink()
}

if (!interactive() && sys.nframe() == 0L &&
    length(commandArgs(trailingOnly = TRUE)) > 0L) main()
