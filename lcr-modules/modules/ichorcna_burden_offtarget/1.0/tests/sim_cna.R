chromosome_bins <- c(
  249, 243, 199, 191, 181, 171, 159, 146, 141, 136, 135,
  134, 115, 107, 102, 90, 83, 80, 59, 64, 47, 51
)

make_genome_skeleton <- function(drop_fraction = 0.12, seed = 1L) {
  chr <- rep(as.character(seq_along(chromosome_bins)), chromosome_bins)
  position <- unlist(lapply(chromosome_bins, function(count) seq_len(count)),
                     use.names = FALSE)
  skeleton <- data.frame(
    chr = chr,
    start = (position - 1) * 1e6,
    end = position * 1e6,
    stringsAsFactors = FALSE
  )
  if (drop_fraction > 0) {
    set.seed(seed)
    keep <- runif(nrow(skeleton)) >= drop_fraction
    skeleton <- skeleton[keep, , drop = FALSE]
  }
  rownames(skeleton) <- NULL
  skeleton
}

expected_logr <- function(f, copy_number) {
  log2((f * copy_number + 2 * (1 - f)) / 2)
}

simulate_cna <- function(f = 0.1, event_count = 12L, sigma = 0.12,
                         seed = 1L, offset = 0, focal = FALSE) {
  genome <- make_genome_skeleton(seed = seed)
  set.seed(seed + 1000L)
  copy_number <- rep(2, nrow(genome))
  chromosomes <- sample(as.character(1:22), event_count, replace = TRUE)
  for (event in seq_len(event_count)) {
    indices <- which(genome$chr == chromosomes[event])
    event_width <- if (focal) {
      max(3L, floor(length(indices) * 0.08))
    } else {
      max(10L, floor(length(indices) * 0.45))
    }
    start <- sample(seq_len(length(indices) - event_width + 1L), 1L)
    affected <- indices[start:(start + event_width - 1L)]
    copy_number[affected] <- sample(c(1, 3, 4), 1L)
  }
  truth <- expected_logr(f, copy_number)
  truth <- truth - median(truth)
  observed <- truth + rnorm(length(truth), sd = sigma) + offset
  data.frame(genome, logR = observed, truth = truth,
             copy_number = copy_number, stringsAsFactors = FALSE)
}

simulate_flat <- function(sigma = 0.12, seed = 1L, offset = 0) {
  simulate_cna(f = 0, event_count = 0L, sigma = sigma,
               seed = seed, offset = offset)
}

truth_amplitude <- function(simulation) {
  sd(simulation$truth - median(simulation$truth))
}

write_cna_seg <- function(simulation, path, sample = "sample") {
  table <- data.frame(
    chr = simulation$chr,
    start = simulation$start,
    end = simulation$end,
    value = simulation$logR,
    stringsAsFactors = FALSE
  )
  names(table)[4L] <- paste0(sample, ".logR")
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

write_corrected_depth <- function(simulation, path) {
  table <- data.frame(
    chr = simulation$chr,
    start = simulation$start,
    end = simulation$end,
    log2_TNratio_corrected = simulation$logR,
    stringsAsFactors = FALSE
  )
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

write_seg_txt <- function(simulation, path, sample = "sample") {
  runs <- rle(simulation$copy_number)
  last <- cumsum(runs$lengths)
  first <- c(1L, head(last, -1L) + 1L)
  table <- data.frame(
    ID = sample,
    chrom = simulation$chr[first],
    start = simulation$start[first],
    end = simulation$end[last],
    num.mark = runs$lengths,
    seg.mean = vapply(seq_along(first), function(index) {
      median(simulation$logR[first[index]:last[index]])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}
