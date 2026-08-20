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

labelled_value <- function(lines, label) {
  prefix <- paste0(label, ":")
  matches <- which(startsWith(lines, prefix))
  if (length(matches) == 0L) {
    return(NA_character_)
  }
  trimws(substring(lines[matches[1L]], nchar(prefix) + 1L))
}

numeric_value <- function(value) {
  suppressWarnings(as.numeric(value))
}

path_metadata <- function(path) {
  components <- strsplit(
    normalizePath(path, mustWork = FALSE), .Platform$file.sep, fixed = TRUE
  )[[1L]]
  seq_build <- components[grepl("^[^-]+--[^-]+$", components)]
  bin_directory <- components[grepl("^bin[0-9.]+kb$", components)]
  seq_build <- if (length(seq_build) > 0L) tail(seq_build, 1L) else NA_character_
  bin_directory <- if (length(bin_directory) > 0L) {
    tail(bin_directory, 1L)
  } else {
    NA_character_
  }
  data.frame(
    seq_type = if (is.na(seq_build)) NA_character_ else sub("--.*$", "", seq_build),
    genome_build = if (is.na(seq_build)) {
      NA_character_
    } else {
      sub("^.*--", "", seq_build)
    },
    binsize_kb = if (is.na(bin_directory)) {
      NA_real_
    } else {
      numeric_value(sub("kb$", "", sub("^bin", "", bin_directory)))
    },
    stringsAsFactors = FALSE
  )
}

read_ichorcna_params <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- which(trimws(lines) == "Sample\tTumor Fraction\tPloidy")
  if (length(header) == 0L || header[1L] == length(lines)) {
    stop("Could not find selected-solution table in params file: ", path,
         call. = FALSE)
  }
  selected <- strsplit(lines[header[1L] + 1L], "\t", fixed = TRUE)[[1L]]
  if (length(selected) < 3L || !nzchar(selected[1L])) {
    stop("Malformed selected-solution row in params file: ", path,
         call. = FALSE)
  }
  metadata <- path_metadata(path)
  cbind(
    data.frame(
      sample_id = selected[1L],
      tumor_fraction = numeric_value(selected[2L]),
      ploidy = numeric_value(selected[3L]),
      gender = labelled_value(lines, "Gender"),
      subclone_fraction = numeric_value(labelled_value(
        lines, "Subclone Fraction"
      )),
      fraction_genome_subclonal = numeric_value(labelled_value(
        lines, "Fraction Genome Subclonal"
      )),
      fraction_cna_subclonal = numeric_value(labelled_value(
        lines, "Fraction CNA Subclonal"
      )),
      coverage = numeric_value(labelled_value(lines, "Coverage")),
      chry_coverage_fraction = numeric_value(labelled_value(
        lines, "ChrY coverage fraction"
      )),
      gc_map_correction_mad = numeric_value(labelled_value(
        lines, "GC-Map correction MAD"
      )),
      stringsAsFactors = FALSE
    ),
    metadata,
    data.frame(
      param_file = normalizePath(path, mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  )
}

aggregate_params <- function(input_dir, output_path) {
  paths <- sort(list.files(
    input_dir, pattern = "[.]param[.]txt$", recursive = TRUE,
    full.names = TRUE
  ))
  if (length(paths) == 0L) {
    stop("No .param.txt files found under: ", input_dir, call. = FALSE)
  }
  table <- do.call(rbind, lapply(paths, read_ichorcna_params))
  duplicated_ids <- unique(table$sample_id[duplicated(table$sample_id)])
  if (length(duplicated_ids) > 0L) {
    stop("Duplicate sample_id values in params files: ",
         paste(duplicated_ids, collapse = ", "), call. = FALSE)
  }
  table <- table[order(table$sample_id), , drop = FALSE]
  rownames(table) <- NULL
  output_directory <- dirname(output_path)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  write.table(table, output_path, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  invisible(table)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_aggregate_args(args)
  required <- c("input_dir", "output")
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
  aggregate_params(options$input_dir, options$output)
}

if (!interactive() && sys.nframe() == 0L &&
    length(commandArgs(trailingOnly = TRUE)) > 0L) main()
