read_burden_tables <- function(paths) {
  if (length(paths) == 0L) {
    stop("No burden TSV files were provided", call. = FALSE)
  }
  tables <- lapply(paths, function(path) {
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
               comment.char = "", quote = "", stringsAsFactors = FALSE)
  })
  shared <- Reduce(intersect, lapply(tables, names))
  if (length(shared) == 0L) {
    stop("Input burden TSVs have no columns in common", call. = FALSE)
  }
  if (any(vapply(tables, function(table) !identical(names(table), shared),
                 logical(1)))) {
    warning("Input TSV schemas differ; taking the column intersection",
            call. = FALSE)
  }
  do.call(rbind, lapply(tables, function(table) table[, shared, drop = FALSE]))
}

rank01 <- function(x) {
  output <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  if (!any(valid)) {
    return(output)
  }
  if (sum(valid) == 1L) {
    output[valid] <- 0.5
  } else {
    output[valid] <- (rank(x[valid], ties.method = "average") - 1) /
      (sum(valid) - 1)
  }
  output
}

add_cohort_annotations <- function(table) {
  required <- c("burden_amp", "sigma_noise", "detect_p")
  missing <- setdiff(required, names(table))
  if (length(missing) > 0L) {
    stop("Missing required burden columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  table$burden_amp_rank <- rank01(as.numeric(table$burden_amp))
  table$sigma_noise_rank <- rank01(as.numeric(table$sigma_noise))
  valid_noise <- as.numeric(table$sigma_noise)
  valid_noise <- valid_noise[is.finite(valid_noise)]
  table$cohort_sigma_noise_p50 <- if (length(valid_noise)) {
    as.numeric(quantile(valid_noise, 0.50, names = FALSE))
  } else {
    NA_real_
  }
  table$cohort_sigma_noise_p90 <- if (length(valid_noise)) {
    as.numeric(quantile(valid_noise, 0.90, names = FALSE))
  } else {
    NA_real_
  }
  p_values <- as.numeric(table$detect_p)
  table$detect_q <- NA_real_
  valid_p <- is.finite(p_values)
  table$detect_q[valid_p] <- p.adjust(p_values[valid_p], method = "BH")
  table$detectable_fdr10 <- is.finite(table$detect_q) & table$detect_q < 0.10
  table
}

qc_summary <- function(table) {
  numeric_names <- names(table)[vapply(table, is.numeric, logical(1))]
  quantile_rows <- lapply(numeric_names, function(name) {
    values <- table[[name]][is.finite(table[[name]])]
    probabilities <- c(0, 0.25, 0.5, 0.75, 0.9, 1)
    quantiles <- if (length(values)) {
      as.numeric(quantile(values, probabilities, names = FALSE))
    } else {
      rep(NA_real_, length(probabilities))
    }
    data.frame(
      summary_type = "statistic",
      metric = name,
      n = length(values),
      q0 = quantiles[1L],
      q25 = quantiles[2L],
      q50 = quantiles[3L],
      q75 = quantiles[4L],
      q90 = quantiles[5L],
      q100 = quantiles[6L],
      value = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  output <- if (length(quantile_rows)) do.call(rbind, quantile_rows) else NULL

  if ("qc_flag" %in% names(table)) {
    flags <- unlist(strsplit(as.character(table$qc_flag), ";", fixed = TRUE))
    counts <- sort(table(flags), decreasing = TRUE)
    flag_rows <- data.frame(
      summary_type = "qc_flag",
      metric = names(counts),
      n = as.integer(counts),
      q0 = NA_real_, q25 = NA_real_, q50 = NA_real_, q75 = NA_real_,
      q90 = NA_real_, q100 = NA_real_, value = NA_real_,
      stringsAsFactors = FALSE
    )
    output <- rbind(output, flag_rows)
  }
  if ("detectable" %in% names(table)) {
    detectable <- as.logical(table$detectable)
    rate <- mean(detectable, na.rm = TRUE)
    if (is.nan(rate)) rate <- NA_real_
    output <- rbind(output, data.frame(
      summary_type = "detection_rate", metric = "detectable", n = sum(!is.na(detectable)),
      q0 = NA_real_, q25 = NA_real_, q50 = NA_real_, q75 = NA_real_,
      q90 = NA_real_, q100 = NA_real_, value = rate,
      stringsAsFactors = FALSE
    ))
  }
  rownames(output) <- NULL
  output
}

aggregate_burden <- function(paths, output_path = NULL, qc_path = NULL) {
  table <- add_cohort_annotations(read_burden_tables(paths))
  summary <- qc_summary(table)
  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    write.table(table, output_path, sep = "\t", quote = FALSE,
                row.names = FALSE, na = "NA")
  }
  if (!is.null(qc_path)) {
    dir.create(dirname(qc_path), recursive = TRUE, showWarnings = FALSE)
    write.table(summary, qc_path, sep = "\t", quote = FALSE,
                row.names = FALSE, na = "NA")
  }
  list(table = table, qc = summary)
}

parse_aggregate_args <- function(args) {
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
      parsed[[pieces[1L]]] <- paste(pieces[-1L], collapse = "=")
    } else {
      if (index == length(args) || startsWith(args[index + 1L], "--")) {
        stop("Missing value for --", token, call. = FALSE)
      }
      index <- index + 1L
      parsed[[token]] <- args[index]
    }
    index <- index + 1L
  }
  parsed
}

aggregate_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_aggregate_args(args)
  required <- c("input_dir", "output", "qc_output")
  missing <- required[!required %in% names(options)]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(paste0("--", missing), collapse = ", "), call. = FALSE)
  }
  paths <- list.files(options$input_dir, pattern = "\\.tsv$", recursive = TRUE,
                      full.names = TRUE)
  paths <- setdiff(normalizePath(paths, mustWork = FALSE),
                   normalizePath(c(options$output, options$qc_output),
                                 mustWork = FALSE))
  aggregate_burden(paths, options$output, options$qc_output)
  invisible(NULL)
}

if (!interactive() && sys.nframe() == 0L &&
    length(commandArgs(trailingOnly = TRUE)) > 0L) aggregate_main()
