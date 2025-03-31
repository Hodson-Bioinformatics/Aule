suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
snakemake@source("PV_utilis.R")

### HELPER FUNCTIONS

get_panel_positions <- function(panel, BSgenome, flank = 300){
  
  # Load panel 
  panel_gr <- import(panel, format = "BED")
  panel_flanked <- panel_gr + flank 
  
  # Extract individual positions
  panel_flanked <- panel_flanked %>% 
    as.data.frame() %>%
    rowwise() %>%
    mutate(start = list(seq(start, end))) %>%
    unnest(start) %>%
    dplyr::select(seqnames, start) %>%
    mutate(end = start)
  
  panel_flanked <- panel_flanked %>%
    mutate(ref_pv1 = as.character(Biostrings::getSeq(BSgenome, 
                                                     GRanges(panel_flanked)))) 
  return(panel_flanked)
}

# Calculate mean PVAF per genotype
calculate_mean_PVAF <- function(tracked_tab){
  
  if(nrow(tracked_tab) == 0){
    data.frame(sample_id = c(),
               tracked_genotype = c(),
               mean_PVAF = c(),
               metric = c())
  } else {
    tracked_tab %>%
      
      dplyr::select(tracked_genotype, contains("PVAF_")) %>%
      gather(sample_id, PVAF, -tracked_genotype) %>%
      mutate(PVAF = as.numeric(PVAF)) %>%
      group_by(sample_id, tracked_genotype) %>%
      summarise(mean_PVAF = mean(PVAF, na.rm = T)) %>%
      mutate(metric = ifelse(grepl("PVAF_total", sample_id), "PVAF_total", "PVAF_qual"),
             sample_id = case_when(sample_id == "PVAF_total" | sample_id == "PVAF_qual" ~ tracked_genotype,
                                   TRUE ~ splitvec(sample_id, "[.]", 1))) %>%
      ungroup() %>%
      distinct()
  }
  
}

# Compute distance between PVs, consider doublets (pv1, pv2) and triplets (pv1, pv2, pv3)
find_alternative_pvs <- function(regions, pv_tab, pv_type){
  pvs_alt <- regions %>%
    dplyr::filter(seqnames == pv_tab$seqnames) %>%
    mutate(ref_pv2 = lead(ref_pv1, pv_tab$dist),
           start_pv2 = lead(start, pv_tab$dist)) %>%
    dplyr::rename(start_pv1 = start) %>%
    dplyr::filter(ref_pv1 == toupper(pv_tab$ref_pv1) & 
                  ref_pv2 == toupper(pv_tab$ref_pv2)) %>%
    dplyr::select(seqnames, start_pv1, start_pv2, ref_pv1, ref_pv2) %>%
    mutate(PV_ID = paste(seqnames, start_pv1, start_pv2, pv_tab$alt_pv1, pv_tab$alt_pv2, sep = "_"))
  
  if (pv_type == "triplets"){
    pvs_alt <- pvs_alt %>%
      mutate(ref_pv3 = lead(ref_pv1, pv_tab$dist2),
             start_pv3 = lead(start_pv1, pv_tab$dist2)) %>%
      dplyr::filter(ref_pv3 == toupper(pv_tab$ref_pv3)) %>%
      dplyr::select(seqnames, start_pv1, start_pv2, start_pv3, ref_pv1, ref_pv2, ref_pv3)  %>%
      mutate(PV_ID = paste(seqnames, start_pv1, start_pv2, start_pv3, pv_tab$alt_pv1, pv_tab$alt_pv2, pv_tab$alt_pv3, sep = "_"))
  }
  return(pvs_alt)
}

# Select alternative PVs
select_alternative_pvs <- function(altenatives){
  if (nrow(altenatives) > 50){
    altenatives[sample(1:nrow(altenatives), size = 50, replace = F),]
  } else {
    size <- (nrow(altenatives)*0.75)+1
    altenatives[sample(1:nrow(altenatives), size = size, replace = T),]
  }
} 

# Find AD
find_supporting_reads <- function(PV_IDs, all_pv){
  ADs <- rep(0, times = length(PV_IDs))
  matches <- PV_IDs %in% all_pv$PV_ID
  if (sum(matches) == 0){
    return(ADs)
  } else {
    ids <- which(matches)
    ADs[ids] <- PV_IDs[ids] %>% 
      map_dbl(~all_pv$AD[.x == all_pv$PV_ID])
    return(ADs)
  }
}

# Compute DP
compute_DP <- function(bam, gr, min_mq = 20){
  param_bam <- ScanBamParam(what = c("qname"),
                            mapqFilter = min_mq,
                            flag = scanBamFlag(isPaired = T,
                                               isProperPair = T,
                                               isDuplicate = F),
                            which = gr)
  reads <- scanBam(bam, param = param_bam, asMates=TRUE) # This is essential to ensure that only unique fragments (read pairs) are counted
  return(length(unique(unlist(reads))))
}

# Compute AF
compute_AF <- function(PV_IDs, ADs, bam, seqnames, end5, end3, min_mapping_qual){
  AFs <- rep(0, times = length(PV_IDs))
  non_zero <- ADs > 0
  if (sum(non_zero) == 0){
    return(AFs)
  } else {
    ids <- which(non_zero)
    AFs[ids] <- ids %>% 
      map_dbl(~compute_DP(bam, GRanges(seqnames = seqnames[.x], 
                                       IRanges(start = end5[.x],
                                               end = end3[.x])),
                          min_mapping_qual))
    AFs[ids] <- ADs[ids]/AFs[ids]
    return(AFs)
  }
}

# Function for selecting a random alternative PV
select_random_PV <- function(pv_tab){
  if (nrow(pv_tab) == 0){
    af = 0
  } else {
    selected <- pv_tab[sample(1:nrow(pv_tab), size = 1),]
    af <- as.numeric(selected$AF)
  }
  return(af)
}

# Compute mean alternative AF fo Monte Carlo interation
mean_alt_af <- function(i, alt_pv_annotated){
  mean_alt_af <- 1:length(alt_pv_annotated) %>%
    map_dbl(~select_random_PV(alt_pv_annotated[[.x]])) %>%
    mean()
}

### MAIN WORKFLOW

mrd_pv <- function(pv_tracked_in, pv_all, pv_type, panel, bams, monte_carlo_iterations, min_mapping_qual){
  
  message(Sys.time(), " Started")
  
  # Load tracked PVs
  load(pv_tracked_in) # Loads tracked_pv
  
  # Calculate mean PVAF per timepoint 
  save(tracked_pv, file = "tracked_pv.RData")
  
  PVAF_tab <- tracked_pv %>%
    map(calculate_mean_PVAF) 
  
  MRD_tab <- PVAF_tab %>%
    map(~.x %>%
          as.data.frame() %>%
          mutate(p_value = NA))
  
  print(MRD_tab)
  PV_alts_tab <- list()
  
  # Define panel positions
  panel_flanked <- get_panel_positions(panel, BSgenome.Hsapiens.UCSC.hg38, flank = 300)
  
  # Link to pv_all 
  names(pv_all) <- basename(pv_all) %>%
    map_chr(~strsplit(.x, split = "[.]")[[1]][1])
  
  # Link to bams 
  names(bams) <- basename(bams) %>%
    map_chr(~strsplit(.x, split = "[.]")[[1]][1])
  
  # Define genotypes and samples
  genotypes <- names(tracked_pv)
  samples <- PVAF_tab %>%
    map(~as.character(.x$sample_id)) %>%
    unlist() %>%
    unique()
  
  message(Sys.time(), " Finished loading")
  # Iterate over genotypes 
  for (genotype in genotypes){
    
    message(Sys.time(), " Processing ", genotype, " genotype ")
    
    # Pre-assign results table
    mrd_genotype <- MRD_tab[[genotype]]
    
    # Find PVs genotype
    genotyped_pv <- tracked_pv[[genotype]]
    
    if (nrow(genotyped_pv) == 0){
      print(genotype)
      MRD_tab[[genotype]] <- data.frame()
      PV_alts_tab[[genotype]] <- data.frame()
      
    } else {
    
    # Define distances between PVs
    genotyped_pv <- genotyped_pv %>%
      mutate(dist = start_pv2 - start_pv1)
    if (pv_type == "triplets"){
      genotyped_pv <- genotyped_pv %>%
        mutate(dist2 = start_pv3 - start_pv1)
    }
    
    PV_alts_tab[[genotype]] <- genotyped_pv %>%
      dplyr::select(PV_ID)
    
     # Find alternative alleles for each PV, save the alternative pvs in a list 
    alt_pv <- 1:nrow(genotyped_pv) %>%
      future_map(~find_alternative_pvs(panel_flanked, genotyped_pv[.x,], pv_type), seed = T) %>%
      # Select 50 random alternative PVs, name them 
      map(~select_alternative_pvs(.x)) 
    names(alt_pv) <- genotyped_pv$PV_ID
    message(Sys.time(), " Alternative PVs found  ")
    
    # Annotate alternative PVs with AD amd DP per timepoint (from pv_all)
    for (sample in samples){
      
      message(Sys.time(), " Processing  ", sample, " sample")
      # Find read pairs containing the alternative PVs
      if (pv_type == "doublets"){
        all_pv <- suppressMessages(read_delim(pv_all[[sample]], delim = "\t")) %>%
          mutate(PV_ID = paste(seqnames, start_pv1+1, start_pv2+1, alt_pv1, alt_pv2, sep = "_")) %>%
          group_by(PV_ID) %>%
          summarise(AD = length(unique(rid)))
      } else if (pv_type == "triplets"){
        all_pv <- suppressMessages(read_delim(pv_all[[sample]], delim = "\t")) %>%
          mutate(PV_ID = paste(seqnames, start_pv1+1, start_pv2+1, start_pv3+1, alt_pv1, alt_pv2, alt_pv3, sep = "_")) %>%
          group_by(PV_ID) %>%
          summarise(AD = length(unique(rid)))
      }
      
      alt_pv_annotated <- alt_pv %>%
        future_map(~.x %>% 
                     mutate(AD = find_supporting_reads(PV_ID, all_pv)))
      message(Sys.time(), " Read pairs found")
      
      # Find DP for entries with non-zero AD (otherwise, AF is zero anyway)
      bamfile <- bams[sample]
      if (pv_type == "doublets"){
        alt_pv_annotated <- alt_pv_annotated %>%
          future_map(~.x %>% 
                       mutate(AF = compute_AF(PV_ID, AD, bamfile, seqnames, start_pv1, start_pv2, min_mapping_qual))) %>%
          map(~.x %>% 
                dplyr::filter(AF < 0.05))
      } else if (pv_type == "triplets"){
        alt_pv_annotated <- alt_pv_annotated %>%
          future_map(~.x %>% 
                       mutate(AF = compute_AF(PV_ID, AD, bamfile, seqnames, start_pv1, start_pv3, min_mapping_qual))) %>%
          map(~.x %>% 
                dplyr::filter(AF < 0.05))
      }
      message(Sys.time(), " Allele frequency computed")
  
      message(Sys.time(), " Starting Monte Carlo simulation ")
      # Monte Carlo simulation to determine p-value of MRD
      alt_afs <- 1:monte_carlo_iterations %>%
        future_map(~mean_alt_af(.x, alt_pv_annotated = alt_pv_annotated), 
                   .options = furrr_options(seed = TRUE))
      PV_alts_tab[genotype][[sample]] <- alt_afs
      
      message(Sys.time(), " Finished Monte Carlo simulation")
      
      # Compute p-values
      for (row in 1:nrow(mrd_genotype)){
        if (mrd_genotype[row, "sample_id"] == sample){
          pval <- mean(mrd_genotype[row, "mean_PVAF"] <= unlist(alt_afs))
          mrd_genotype[row, "p_value"] <- pval
          message(sample, " p-value: ", pval, " (id) ", Sys.time())
        }
      }
      
      # Update MRD table
      MRD_tab[[genotype]] <- mrd_genotype
      
     }
      
  } 
    }
    
  # Collapse the final table of MRD p-values
  MRD_tab_collapsed <-do.call(rbind, MRD_tab)
  
  message(Sys.time(), " All done! ")
  output <- list(MRD_tab_collapsed = MRD_tab_collapsed,
                 PV_alts_tab = PV_alts_tab)
  
  return(output)
   
}

## Load params
pv_tracked_in <- snakemake@input[["pv_tracked"]]
pv_all        <- snakemake@input[["pv_all"]]
bams          <- snakemake@input[["bams"]]
panel         <- snakemake@input[["panel"]]
pv_type       <- snakemake@params[["pv_type"]]
monte_carlo_iterations <- snakemake@params[["monte_carlo_iterations"]]
min_mapping_qual <- snakemake@params[["min_mapping_qual"]]
pv_mrd        <- snakemake@output[["pv_mrd"]] 
pv_tab        <- snakemake@output[["pv_tab"]] 
threads       <- snakemake@threads

# Register cores for parallel run
plan(multisession, workers = threads)

# # ## Develop
# basedir = "/home/jak75/rds/rds-hodson-rds-elshhQhj4ec/runs/DIRECT/20240605_PV/results/aule_pv-1.0/"
# pv_tracked_in <- paste0(basedir, "07-track/capture--hg38/DIR1195.doublets.tracked_pv.RData")
# file.exists(pv_tracked_in)
# pv_all <- list.files(paste0(basedir, "02-pv_all/capture--hg38/"),
#                      pattern = "^DIR1195.*doublets.all.txt$", full.names = T)
# bams <- list.files(paste0(basedir, "03-sorted_bam/capture--hg38"),
#                    pattern = "^DIR1195.*bam$", full.names = T)
# min_bq <- 20
# min_mq <- 20
# pv_type <- "doublets"
# panel <- "/home/jak75/rds/rds-hodson-rds-elshhQhj4ec/runs/DIRECT/20240605_PV/data/exome_bed/hg38/DIRECT_2021-BED.bed"

# Track PVs
mrd_stat <-  mrd_pv(pv_tracked_in = pv_tracked_in, 
                    pv_all = pv_all, 
                    pv_type = pv_type, 
                    panel = panel, 
                    bams  = bams,
                    monte_carlo_iterations = monte_carlo_iterations,
                    min_mapping_qual = min_mapping_qual)

mrd_stat_tab <- mrd_stat$MRD_tab_collapsed
mrd_pv_tab <- mrd_stat$PV_alts_tab

write_delim(mrd_stat_tab, file = pv_mrd, delim = "\t")
save(mrd_pv_tab, file = pv_tab)
