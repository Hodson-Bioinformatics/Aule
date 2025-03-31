##### Script ot perform basic filtering (R1/R2 concordance, MQ, BQ)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))
snakemake@source("PV_utilis.R")

load_and_collapse <- function(pv_files){
  n_total <- length(pv_files)
  pv_long <- pv_files
  pv_long <- pv_long  %>%
    future_map(~ suppressMessages(read_delim(.x, delim = "\t")) %>%
                 dplyr::select(PV_ID, AD_total, contains("DP_")) %>%
                 rowwise() %>%
                 mutate(mean_DP = mean(c_across(contains("DP_"))),
                        PVAF = AD_total/mean_DP)) %>%
    bind_rows() %>%
    dplyr::filter(!is.infinite(PVAF)) %>%
    mutate(PVAF = ifelse(PVAF > 1, 1, PVAF)) %>%
    group_by(PV_ID) %>%
    summarise(AD_total = sum(AD_total, na.rm = T),
              DP_total = round(sum(mean_DP, na.rm = T), 2), 
              AD_average = mean(AD_total, na.rm = T),
              N = length(PV_ID),
              maxPVAF = max(PVAF, na.rm = T),
              meanPVAF = mean(PVAF, na.rm = T)) %>%
    mutate(PVAF_total = AD_total/DP_total, 
           PVAF_total = ifelse(PVAF_total > 1, 1, PVAF_total)) %>%
    arrange(desc(N), desc(PVAF_total))
  return(pv_long)
}


## Load params
doublets    <- snakemake@input[["doublets"]]
triplets    <- snakemake@input[["triplets"]]
pon_samples <- snakemake@params[["pon_samples"]]
panel       <- snakemake@input[["panel"]]
threads     <- snakemake@threads
output_doublets      <- snakemake@output[["pon_doublets"]]
output_triplets     <- snakemake@output[["pon_triplets"]]

# Load files 
doublets <- list.files("results/aule_pv_pon-1.0/04-pv_annotated/capture--hg38", pattern = "*doublets.annotated.txt", full.names = T)
triplets <- list.files("results/aule_pv_pon-1.0/04-pv_annotated/capture--hg38", pattern = "*triplets.annotated.txt", full.names = T)
names(doublets) <- gsub(".doublets.annotated.txt", "", basename(doublets))
names(triplets) <- gsub(".triplets.annotated.txt", "", basename(triplets))

# Filter by the samples list 
if (!file.exists(pon_samples)){
  message("No PON samples IDs provided, processing all")
} else {
  pon_samples <- read.delim(pon_samples, col.names = "sample_id")[,1]
  doublets <- doublets[pon_samples]
  triplets <- triplets[pon_samples]
}
# Set multithread run 
plan(multisession, workers = threads)

# Generate PON and save 
message("Summarising doublets")
doublets_summarised <- load_and_collapse(doublets) 
write_delim(doublets_summarised, file = output_doublets, delim = "\t")
rm(doublets_summarised)

message("Summarising triplets")
triplets_summarised <- load_and_collapse(triplets) 
write_delim(triplets_summarised, file = output_triplets, delim = "\t")

