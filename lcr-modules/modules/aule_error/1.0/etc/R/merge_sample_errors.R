# Author: Ilias Moutsopoulos
# Date: 29/01/2024
# Description: Merge errors across samples

## If R is installed with conda, this will configure the 
## library paths to only point to the conda R package directory
local({
  conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
  rhome <- normalizePath(R.home())
  if (!is.na(conda) && startsWith(rhome, conda)) {
    .libPaths(R.home("library"), include.site = FALSE)
  }
})

# Parse snakemake object
errors                 <- snakemake@input[["errors"]]
stats                  <- snakemake@input[["stats"]]
errors_split_by_sample <- snakemake@output[["errors_split_by_sample"]]
errors_merged          <- snakemake@output[["errors_merged"]]
stats_merged           <- snakemake@output[["stats"]]
log_stdout             <- snakemake@log[["stdout"]]
log_stderr             <- snakemake@log[["stderr"]]


# Capture output
file_stderr <- file(log_stderr, open = "at")
sink(file_stderr, type = "message", append = TRUE)
sink(log_stdout, type = "output", append = TRUE)

library(dplyr)

stats_df <- data.frame()
for (i in seq_along(stats)) {
  stats_temp <- read.csv(stats[i]) %>%
    dplyr::mutate(sample = basename(dirname(stats[i])), .before = 1)
  
  stats_df <- rbind(stats_df, stats_temp)
}

coords_concat <- data.frame()
pb = txtProgressBar(min = 0, max = length(errors), initial = 0, style = 3)
for (i in seq_along(errors)) {
  coords_temp <- read.csv(errors[i]) %>%
    dplyr::mutate(sample = basename(dirname(errors[i])), .before = 1)
  
  if (nrow(coords_concat > 0) & nrow(coords_temp) > 0) {
    if (!identical(colnames(coords_concat), colnames(coords_temp))) {
      stop("Input tables do not have identical column names")
    }
  }
  
  coords_concat <- rbind(coords_concat, coords_temp)
  
  setTxtProgressBar(pb,i)
}

coords_merged <- coords_concat %>%
  dplyr::group_by(index, name, unique_name, chr, pos, chr_pos, seq1, seq3, seq5, seq7) %>%
  dplyr::summarise(
    reads                     = sum(reads),
    families                  = sum(families),
    reads_nontrivial          = sum(reads_nontrivial),
    families_nontrivial       = sum(families_nontrivial),
    families_size1            = sum(families_size1),
    families_size2            = sum(families_size2),
    families_size3            = sum(families_size3),
    r12_overlap_depth         = sum(r12_overlap_depth),
    total_base_quality        = sum(total_base_quality),
    total_mapping_quality     = sum(total_mapping_quality),
    total_distance_from_5p    = sum(total_distance_from_5p),
    total_distance_from_3p    = sum(total_distance_from_3p),
    total_fragment_length     = sum(total_fragment_length),
    umi_errors                = sum(umi_errors),
    r12_discordant            = sum(r12_discordant),
    aligned_base_A            = sum(aligned_base_A),
    aligned_base_T            = sum(aligned_base_T),
    aligned_base_C            = sum(aligned_base_C),
    aligned_base_G            = sum(aligned_base_G),
    aligned_base_nontrivial_A = sum(aligned_base_nontrivial_A),
    aligned_base_nontrivial_T = sum(aligned_base_nontrivial_T),
    aligned_base_nontrivial_C = sum(aligned_base_nontrivial_C),
    aligned_base_nontrivial_G = sum(aligned_base_nontrivial_G)
  )

write.csv(coords_concat, errors_split_by_sample, row.names = FALSE)
write.csv(coords_merged, errors_merged, row.names = FALSE)
write.csv(stats_df, stats_merged, row.names = FALSE)

# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)
