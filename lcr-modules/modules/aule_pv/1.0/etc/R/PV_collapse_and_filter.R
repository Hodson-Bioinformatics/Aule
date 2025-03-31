##### Script ot perform basic filtering (R1/R2 concordance, MQ, BQ)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Rsamtools))
snakemake@source("PV_utilis.R")

# Filter PVs by R1/R2, MAPQ and BQ
filter_PV <- function(PVs, min_mq = 0, min_bq = 0, r1r2_mode = "concordant", type = "doublets"){
  
  # # Check if PV is present in R1 and R2 TBA
  # if (r1r2_mode == "strict"){
  #   PVs <- PVs %>%
  #     dplyr::filter(both_pv1 == 1 & both_pv2 == 1)
  #   if (type == "triplets"){
  #     PVs <- PVs  %>%
  #       dplyr::filter(both_pv3 == 1)
  #   }
  # } 
  # For the concordant mode - do not apply separate filtering, R1/R2 concordance has been checked at the PVs finding stage
  
  # Filter by MQ, BQ and family concordance
  PVs <- PVs %>%
    dplyr::filter(mq_r1 >= min_mq & mq_r2 >= min_mq & 
                    concordance_pv1 == 0 & concordance_pv2 == 0) %>%
    mutate(bq_passed = 
             bq_r1_pv1 >= min_bq & bq_r2_pv1 >= min_bq &
             bq_r1_pv2 >= min_bq & bq_r2_pv2 >= min_bq)
  if (type == "triplets"){
    PVs <- PVs %>%
      dplyr::filter(concordance_pv3 == 0) %>%
      mutate(bq_passed = 
               bq_r1_pv1 >= min_bq & bq_r2_pv1 >= min_bq &
               bq_r1_pv2 >= min_bq & bq_r2_pv2 >= min_bq &
               bq_r1_pv3 >= min_bq & bq_r2_pv3 >= min_bq)
  }
  return(PVs)
}

# Collapse PVs by region
collapse_PV <- function(PVs, detailed= F){
  
  # Find columns
  cols <- colnames(PVs)
  start_cols <- cols[grepl("start_", cols)]
  ref_cols <- cols[grepl("ref_", cols)]
  alt_cols <- cols[grepl("alt_", cols)]
  type_cols <- cols[grepl("type_", cols)]
  bq_r1 <- cols[grepl("bq_r1", cols)]
  bq_r2 <- cols[grepl("bq_r2", cols)]
  
  # Summarise
  if(detailed){
    pvs_collapsed <- PVs %>%
      arrange(seqnames, !!!syms(start_cols)) %>%
      group_by(PV_ID, seqnames, !!!syms(start_cols),
               !!!syms(ref_cols), !!!syms(alt_cols),
               !!!syms(type_cols), prime5, prime3, tl, bq_passed, PV_ID_init) %>%
      summarise(AD = length(PV_ID)) %>%
      arrange(desc(AD)) %>%
      mutate(end5 = min(!!!syms(start_cols)),
             end3 = max(!!!syms(start_cols)),
             FILTER = ".")
  } else {
    pvs_collapsed <- PVs %>%
      arrange(seqnames, !!!syms(start_cols)) %>%
      group_by(PV_ID, seqnames, !!!syms(start_cols),
               !!!syms(ref_cols), !!!syms(alt_cols),
               !!!syms(type_cols), PV_ID_init) %>%
      summarise(AD_qual = length(unique(rid[bq_passed])),
                AD_total = length(unique(rid))) %>%
      arrange(desc(AD_qual)) %>%
      mutate(end5 = min(!!!syms(start_cols)),
             end3 = max(!!!syms(start_cols)),
             FILTER = ".")
  }
  return(pvs_collapsed)
}

filter_and_annotate <- function(doublets_all, triplets_all, bam,  min_mq = 0, min_bq = 0, r1r2 = "concordant", detailed = F){

  # Load PV files
  message("Loading PVs")
  doublets <- suppressMessages(read_delim(doublets_all, delim = "\t")) %>%
    mutate(PV_ID_init = paste(seqnames, start_pv1, start_pv2, alt_pv1, alt_pv2, sep = "_"),
           start_pv1 = start_pv1+1,
           start_pv2 = start_pv2+1,
           PV_ID = paste(seqnames, start_pv1, start_pv2, alt_pv1, alt_pv2, sep = "_"))

  triplets <- suppressMessages(read_delim(triplets_all, delim = "\t")) %>%
    mutate(PV_ID_init = paste(seqnames, start_pv1, start_pv2, start_pv3, alt_pv1, alt_pv2, alt_pv3, sep = "_"),
           start_pv1 = start_pv1+1,
           start_pv2 = start_pv2+1,
           start_pv3 = start_pv3+1,
           PV_ID = paste(seqnames, start_pv1, start_pv2, start_pv3, alt_pv1, alt_pv2, alt_pv3, sep = "_"))

  # Filter PVs (present on both reads with satisfactory quality)
  doublets_keep <- filter_PV(doublets, min_mq, min_bq, r1r2, "doublets")
  triplets_keep <- filter_PV(triplets, min_mq, min_bq, r1r2, "triplets")

  # Collapse
  doublets_collapsed <- collapse_PV(doublets_keep, detailed)
  triplets_collapsed <- collapse_PV(triplets_keep, detailed)
  rm(doublets, triplets) # Some housekeeping 
  gc()
  
  # Extract unique genomic positions with PVs 
  pv_ranges <- data.frame(
    seqnames = c(rep(doublets_collapsed$seqnames, 2), rep(triplets_collapsed$seqnames, 3)),
    start = 1, 
    end = c(doublets_collapsed$start_pv1, doublets_collapsed$start_pv2, 
            triplets_collapsed$start_pv1, triplets_collapsed$start_pv2, triplets_collapsed$start_pv3)
  ) %>%
    mutate(end = as.numeric(end), 
           start = end) %>%
    distinct() %>%
    arrange(seqnames, start) %>%
    GRanges()
  
  return(list(doublets = doublets_collapsed,
              triplets = triplets_collapsed,
              ranges = pv_ranges))
}

## Load params
doublets  <- snakemake@input[["doublets"]]
triplets  <- snakemake@input[["triplets"]]
bam       <- snakemake@input[["bam"]]
min_mq    <- snakemake@params[["min_mapping_qual"]]
min_bq    <- snakemake@params[["min_base_qual"]]
r1r2      <- as.logical(snakemake@params[["r1r2"]])
detailed  <- as.logical(snakemake@params[["detailed"]])
doublets_out <- snakemake@output[["doublets"]]
triplets_out <- snakemake@output[["triplets"]]
bed <- snakemake@output[["bed"]]


# Collapse PVs
collapsed <- filter_and_annotate(doublets_all = doublets, 
                                 triplets_all = triplets, 
                                 bam = bam,
                                 min_mq = min_mq, 
                                 min_bq = min_bq, 
                                 r1r2 = r1r2,
                                 detailed = detailed)

write_delim(collapsed$doublets, file = doublets_out, delim = "\t")
write_delim(collapsed$triplets, file = triplets_out, delim = "\t")
export(collapsed$ranges, bed, format = "BED")