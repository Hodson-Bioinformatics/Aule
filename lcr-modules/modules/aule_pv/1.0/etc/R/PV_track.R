##### Script to select PVs for tracking
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
snakemake@source("PV_utilis.R")

track_pv <- function(pv_genotypes, pv_all, pv_type, bams, min_bq, min_mq){
  
  message(Sys.time(), " Started")
  
  # Empty output 
  tracked_pv <- list()
  stats_pv <- list()
  summarised_pv <- list()
    
  # Load genotypes 
  message(Sys.time(), " Load genotypes")
  genotypes_tabs <- pv_genotypes %>%
    map(~read_delim(.x, delim = "\t"))
  names(genotypes_tabs) <- basename(pv_genotypes) %>%
    map_chr(~strsplit(.x, split = "-")[[1]][1])
  
  # Define baselines 
  baseline <- names(genotypes_tabs)
  baseline <- baseline[grepl("Biopsy_B|Plasma_B|Plasma_PD|Biopsy_PD", baseline)]
  
  # Link to pv_all 
  names(pv_all) <- basename(pv_all) %>%
    map_chr(~strsplit(.x, split = "[.]")[[1]][1])
  
  # Link to bams 
  names(bams) <- basename(bams) %>%
    map_chr(~strsplit(.x, split = "[.]")[[1]][1])
  
  # Collect filtering stats
  stats <- genotypes_tabs %>%
    map_dfr(~stat_genotypes(.x), .id = "sample_id") 
  
  # Focus only on PVs passing filters
  message(Sys.time(), " Filter genotypes")
  genotypes_tabs_passed <- genotypes_tabs %>%
    map(~dplyr::filter(.x, FILTER == "PASS"))
  
  # Filter genotypes tabs ()
  # Iterate through the baseline samples 
  for (b in seq_along(baseline)){
  
    baseline_name <- baseline[b]
    baseline_genotype <- genotypes_tabs_passed[[baseline_name]]
    baseline_stats <- list()
    message(Sys.time(), " Processing ", baseline_name, " genotype: ", nrow(baseline_genotype), " PVs found")
    
    print(baseline_genotype[1:5,])
    # Find reads supporting the genotyped PVs in all timepoints 
    timepoints <- names(genotypes_tabs)
    for (t in timepoints){
      message(Sys.time(), " Finding supporting reads in ", t)
      
      # Extract individual metrics of tracked PVs
      genotype_pv <- genotypes_tabs[[t]] 
      if(nrow(baseline_genotype) > 0){
        genotype_pv <- get_genotype_pv(PV_IDs = baseline_genotype$PV_ID, 
                                      genotype_pv = genotype_pv, 
                                      timepoint = t)
      
        # Save it to the tracking table
        baseline_genotype <- baseline_genotype %>%
          left_join(genotype_pv, by = "PV_ID") %>%
          mutate_if(is.numeric, ~tidyr::replace_na(., 0)) 
        
        # Extract number of unique reads supporting the tracked PVs
        if (pv_type == "doublets"){
          all_pv <- read_delim(pv_all[[t]], delim = "\t") %>%
            mutate(PV_ID = paste(seqnames, start_pv1+1, start_pv2+1, alt_pv1, alt_pv2, sep = "_")) %>%
            dplyr::filter(PV_ID %in% genotype_pv$PV_ID) %>%
            mutate(qual_passed = bq_r1_pv1 >= min_bq & bq_r2_pv1 >= min_bq &
                     bq_r1_pv2 >= min_bq & bq_r2_pv2 >= min_bq) %>%
            dplyr::filter(mq_r1 >= min_mq & mq_r2 >= min_mq)
        } else if (pv_type == "triplets"){
          all_pv <- read_delim(pv_all[[t]], delim = "\t") %>%
            mutate(PV_ID = paste(seqnames, start_pv1+1, start_pv2+1, start_pv3+1, alt_pv1, alt_pv2, alt_pv3, sep = "_")) %>%
            dplyr::filter(PV_ID %in% genotype_pv$PV_ID) %>%
            mutate(qual_passed = bq_r1_pv1 >= min_bq & bq_r2_pv1 >= min_bq & 
                     bq_r1_pv2 >= min_bq & bq_r2_pv2 >= min_bq &
                     bq_r1_pv3 >= min_bq & bq_r2_pv3 >= min_bq) %>%
            dplyr::filter(mq_r1 >= min_mq & mq_r2 >= min_mq)
        }
       
        
        # Extract the total number of reads covering the PV region
        all_pv_gr <- pv_ranges(baseline_genotype, combined = T)
        
        param_bam <- ScanBamParam(what = c("qname"),
                                  mapqFilter = min_mq,
                                  flag = scanBamFlag(isPaired = T,
                                                     isProperPair = T,
                                                     isDuplicate = F),
                                  which = all_pv_gr)
        bamfile <- bams[t]
        reads <- scanBam(bamfile, param = param_bam, asMates=TRUE) # This is essential to ensure that only unique fragments (read pairs) are counted
        n_informative_reads <- length(unique(unlist(reads)))
        
        # Summarised stats
        timepoint_stats <- data.frame(
          pv_unique_reads_total = length(unique(all_pv$rid)),
          pv_unique_reads_qual = length(unique(all_pv$rid[all_pv$qual_passed])),
          pv_max_fs = max(c(all_pv$r1_fs, all_pv$r2_fs)),
          pv_nontrivial_fs = length(unique(all_pv$rid[all_pv$r1_fs > 1 |  all_pv$r2_fs > 1])),
          pv_rev = length(unique(all_pv$rid[all_pv$rev_r1 == 1])),
          n_informative_reads = n_informative_reads
        ) %>%
          mutate(
            PVAF_combined_qual = pv_unique_reads_qual/n_informative_reads,
            PVAF_combined_total = pv_unique_reads_total/n_informative_reads,
            sample_id = t
          )
      } else {
        timepoint_stats <- data.frame(
          pv_unique_reads_total = NA,
          pv_unique_reads_qual = NA,
          pv_max_fs = NA,
          pv_nontrivial_fs = NA,
          pv_rev = NA,
          n_informative_reads = NA,
          PVAF_combined_qual = NA,
          PVAF_combined_total = NA,
          sample_id = t)
      }
    
      # Save it to the summary table
      baseline_stats[[t]] <- timepoint_stats
      
    }
    
    # Combine the tracking results
    baseline_genotype <- baseline_genotype %>%
          mutate(tracked_genotype = baseline_name) %>%
          relocate(tracked_genotype) %>%
          relocate(contains("AD_qual_tumour"), 
                   contains("AD_total_tumour"),
                   contains("PVAF_qual"),
                   contains("PVAF_total"),
                   .after = "PVAF_total")
    
    tracked_pv[[baseline_name]] <- baseline_genotype
    
    # Combine the tracking metrics 
    baseline_stats <- baseline_stats %>%
      map_df(~.x)
    summarised_pv[[baseline_name]] <- baseline_stats
      }
  
  # Collapse into a single dataframe
  summarised_pv <- summarised_pv %>%
    map_df(~.x, .id = "patient_id")
  
  return(list(summarised_pv = summarised_pv, 
              tracked_pv = tracked_pv))
}

## Load params
pv_genotypes  <- snakemake@input[["pv_genotypes"]]
pv_all        <- snakemake@input[["pv_all"]]
bams          <- snakemake@input[["bams"]]
pv_type       <- snakemake@params[["pv_type"]]
min_mq        <- snakemake@params[["min_mapping_qual"]]
min_bq        <- snakemake@params[["min_base_qual"]]
pv_tracked    <- snakemake@output[["pv_tracked"]] 
pv_summarised <- snakemake@output[["pv_summarised"]] 

# ## Develop
# pv_genotypes <- list.files("/servers/bio-data-groups-gagarin/hodson/Aule/runs/20231208_test_SD/results/aule_pv-1.0/06-genotype/capture--hg38",
#                            "doublets.genotype.txt", full.names = T)
# pv_genotypes <- pv_genotypes[grepl("DIR1003", pv_genotypes)]
# pv_all <- list.files("/servers/bio-data-groups-gagarin/hodson/Aule/runs/20231208_test_SD/results/aule_pv-1.0/02-pv_all/capture--hg38",
#                            "doublets.all.txt", full.names = T)
# pv_all <- pv_all[grepl("DIR1003", pv_all)]
# pv_collapsed <- list.files("/servers/bio-data-groups-gagarin/hodson/Aule/runs/20231208_test_SD/results/aule_pv-1.0/04-pv_collapsed/capture--hg38",
#                      "doublets.raw.txt", full.names = T)
# pv_collapsed <- pv_collapsed[grepl("DIR1003", pv_collapsed)]
# bams <- list.files("/servers/bio-data-groups-gagarin/hodson/Aule/runs/20231208_test_SD/bams/capture--hg38/",
#                    ".bam$", full.names = T)
# 
# baseline <- c("DIR1003_Biopsy_B_v1", "DIR1003_Plasma_B_v1", "DIR1003_Plasma_PD_v1")
# min_bq <- 20
# min_mq <- 20
# pv_type = "doublets"

# Track PVs
tracked <-  track_pv(pv_genotypes = pv_genotypes, 
                     pv_all = pv_all,
                     pv_type = pv_type,
                     bams = bams,
                     min_bq = min_bq, 
                     min_mq = min_mq)

# Save results
summarised_pv <- tracked$summarised_pv
tracked_pv <- tracked$tracked_pv
write_delim(summarised_pv, file = pv_summarised, delim = "\t")
save(tracked_pv, file = pv_tracked)


