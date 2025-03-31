### HELPER FUNCTIONS
# Helper functions for PV analysis

# # Vectrorized strsplit
splitvec <- function(vector, split, select, merge = "_"){
  processed <- sapply(vector, function(x){
    separated <- unlist(strsplit(x, split = split))[select]
    if (length(separated) > 1){
      return(paste(separated, collapse = merge))
    } else
      return(separated)
  })
  processed <- unname(processed)
  return(processed)
}

# PVs as GRanges (in 1-based coordinates)
granges_PV <- function(PV_file = doublets, type = "doublets"){
  if(type == "doublets"){
    gr <- GRanges(seqnames = PV_file$seqnames,
                  ranges = IRanges(start = PV_file$start_pv1,
                                   end = PV_file$start_pv2)) %>%
      reduce()
  } else if (type == "triplets"){
    PV_file <- PV_file %>%
      mutate(end5 = min(start_pv1, start_pv2, start_pv3),
             end3 = max(start_pv1, start_pv2, start_pv3))
    gr <- GRanges(seqnames = PV_file$seqnames,
                  ranges = IRanges(start = PV_file$end5,
                                   end = PV_file$end3)) %>%
      reduce()
  }
  return(gr)
}

# Function to return loci participating in PVs as GRanges
loci_PV <- function(doublets, triplets){
  phased_loci <- c(GRanges(seqnames = doublets$seqnames,
                           ranges = IRanges(start = doublets$start_pv1, width = 1)),
                   GRanges(seqnames = doublets$seqnames,
                           ranges = IRanges(start = doublets$start_pv2, width = 1)),
                   GRanges(seqnames = triplets$seqnames,
                           ranges = IRanges(start = triplets$start_pv1, width = 1)),
                   GRanges(seqnames = triplets$seqnames,
                           ranges = IRanges(start = triplets$start_pv2, width = 1)),
                   GRanges(seqnames = triplets$seqnames,
                           ranges = IRanges(start = triplets$start_pv3, width = 1))) %>%
    reduce()
  return(phased_loci)
}

# Function to return basic metrics for PVs
metrics_PV <- function(doublets, triplets, bam){

  #Compute the total number of PVS
  uq_doublets <- length(unique(doublets$PV_ID))
  uq_triplets <- length(unique(triplets$PV_ID))

  n_doublets <- length(doublets$PV_ID)
  n_triplets <- length(triplets$PV_ID)

  #Compute the total number of PVS
  uqR_doublets <- length(unique(doublets$rid))
  uqR_triplets <- length(unique(triplets$rid))

  #Compose GRanges objects
  doublets_gr <- granges_PV(doublets, "doublets")
  triplets_gr <- granges_PV(triplets, "triplets")

  # Compute PV loci
  totalBases <- sum(width(bam))
  phased_loci <- loci_PV(doublets, triplets)

  #Compute the number of IR
  ir_doublets <- length(unique(subjectHits(findOverlaps(doublets_gr, bam))))
  ir_triplets <- length(unique(subjectHits(findOverlaps(triplets_gr, bam))))

  return(list(doublets = list(ir = ir_doublets,
                              uq = uq_doublets,
                              n = n_doublets,
                              uqR = uqR_doublets),
              triplets = list(ir = ir_triplets,
                              uq = uq_triplets,
                              n = n_triplets,
                              uqR = uqR_triplets),
              totalBases = totalBases,
              #phased_loci = phased_loci,
              n_phased_loci = length(phased_loci)))
}

# Load panel
load_panel <- function(panel){
  panel <- suppressMessages(read_delim(panel, delim = "\t",
                      col_names = c("seqnames", "start", "end", "name", "ID"))) %>%
  group_by(seqnames, start, end) %>%
  summarise(name = name[1],
            start = start+1,
            ID = ID[1])
  panel_gr <- GRanges(panel)
  return(panel_gr)
}

# Generate PVs genomic ranges
pv_ranges <- function(PVs, combined = F){
  
  if (combined){
    gr_tab <- data.frame(seqnames = rep(PVs$seqnames, 2),
                         start = c(PVs$start_pv1, PVs$start_pv2)-1,
                         end = c(PVs$start_pv1, PVs$start_pv2),
                         ALT = c(PVs$alt_pv1, PVs$alt_pv2))
    if("start_pv3" %in%  colnames(PVs)){
      gr_tab <- rbind(gr_tab,
                      data.frame(seqnames = PVs$seqnames,
                                 start = PVs$start_pv3-1,
                                 end = PVs$start_pv3,
                                 ALT = PVs$alt_pv3))
    }
    gr_tab <- gr_tab %>% distinct()
    gr_tab <- GRanges(gr_tab)
    return(gr_tab)
    
  } else {
    pv1 <- GRanges(seqnames = PVs$seqnames,
                   ranges = IRanges(start = PVs$start_pv1, width = 1),
                   ALT = PVs$alt_pv1)
    pv2 <- GRanges(seqnames = PVs$seqnames,
                   ranges = IRanges(start = PVs$start_pv2, width = 1),
                   ALT = PVs$alt_pv2)
    if("start_pv3" %in%  colnames(PVs)) {
      pv3 <- GRanges(seqnames = PVs$seqnames,
                     ranges = IRanges(start = PVs$start_pv3, width = 1),
                     ALT = PVs$alt_pv3)
      return(list(pv1 = pv1, pv2 = pv2, pv3 = pv3))
    } else {
      return(list(pv1 = pv1, pv2 = pv2))
    }
  }
}

#### GENOTYPING
flag_germline <- function(PVs, germline){
  # Load germline PVs
  germline <- germline %>%
    mutate(Sample_type = "Germline")

  PVs <- PVs %>%
    mutate(Germline = PV_ID %in% germline$PV_ID,
           FILTER = case_when(Germline & FILTER == "." ~ "germline",
                              Germline ~ paste(FILTER, "germline", sep = ";"),
                              TRUE ~ FILTER))
  return(PVs)
}

flag_too_low <- function(PVs, minDP, minAF, minAD){

  PVs %>%
    mutate(FILTER = case_when(DP_pv < minDP & FILTER == "." ~ "too_low_DP",
                              DP_pv < minDP ~ paste(FILTER, "too_low_DP", sep = ";"),
                              TRUE ~ FILTER),
           FILTER = case_when(AD < minAD & FILTER == "." ~ "too_low_AD",
                              AD < minAD ~ paste(FILTER, "too_low_AD", sep = ";"),
                              TRUE ~ FILTER),
           FILTER = case_when(AF_pv < minAF & FILTER == "." ~ "too_low_AF",
                              AF_pv < minAF ~ paste(FILTER, "too_low_AF", sep = ";"),
                              TRUE ~ FILTER))
}

find_tumour_markers <- function(ffpe, plasma, ID, Unit){
 rbind(ffpe, plasma) %>%
    dplyr::filter(FILTER == ".") %>%
    mutate(FILTER = "PASS",
           ID = ID,
           Unit = Unit) %>%
    arrange(PV_ID) %>%
    group_by(PV_ID) %>%
    mutate(Genotyping = case_when(length(PV_ID) == 2 ~ "FFPE_and_plasma",
                                  Sample_type == "Baseline_tumour" ~ "FFPE_only",
                                  Sample_type == "Baseline_plasma" ~ "Plasma_only"))

}

###### MONITORING
getPVcoverage <- function(bamfile, region, minMQ){
  param <- ScanBamParam(which=region,
                        what = c("qname", "mapq"),
                        mapqFilter = 0,
                        flag = scanBamFlag(isPaired = T,
                                           isProperPair = T,
                                           isDuplicate = F))
  reads <- scanBam(bamfile, param = param, asMates=TRUE) # This is essential to ensure that only unique fragments (read pairs) are counted
  reads <- lapply(reads, function(x){data.frame(qname = x$qname,
                                                mapq = x$mapq)})
  reads <- do.call("rbind", reads)
  bam_stats <- tibble(all_IR = length(unique(reads$qname)),
                      qual_IR = length(unique(reads$qname[reads$mapq >= minMQ])))
  return(bam_stats)
}

getPVreads <- function(PVs_all_file, PVs, minMQ, minBQ, PV_type){

  if(PV_type == "doublets"){
    PVs_all <- suppressMessages(read_delim(PVs_all_file, delim = "\t")) %>%
      mutate(start_pv1 = start_pv1+1,
             start_pv2 = start_pv2+1,
             PV_ID = paste(seqnames, start_pv1, start_pv2, alt_pv1, alt_pv2, sep = "_")) %>%
      dplyr::filter(PV_ID %in% PVs$PV_ID & both_pv1 == 1 & both_pv2 == 1 &
                      bq_r1_pv1 > minBQ & bq_r1_pv2 > minBQ & bq_r2_pv1 > minBQ & bq_r2_pv2 > minBQ) %>%
      arrange(PV_ID)
  } else {
    PVs_all <- suppressMessages(read_delim(PVs_all_file, delim = "\t")) %>%
      mutate(start_pv1 = start_pv1+1,
             start_pv2 = start_pv2+1,
             start_pv3 = start_pv3+1,
             PV_ID = paste(seqnames, start_pv1, start_pv2, start_pv3, alt_pv1, alt_pv2, alt_pv3, sep = "_")) %>%
      dplyr::filter(PV_ID %in% PVs$PV_ID & both_pv1 == 1 & both_pv2 == 1 &
                      bq_r1_pv1 > minBQ & bq_r1_pv2 > minBQ & bq_r1_pv3 > minBQ & bq_r2_pv1 > minBQ & bq_r2_pv2 > minBQ & bq_r2_pv3 > minBQ) %>%
      arrange(PV_ID)
  }

  read_stats <- tibble(total_PV = length(unique(PVs_all$rid)),
                       qual_PV = length(unique(PVs_all$rid[PVs_all$mq_r1 >= minMQ |
                                                             PVs_all$mq_r2 >= minMQ])))
  return(read_stats)
}

# Function to extract filtering stats from genotype table 
stat_genotypes <- function(pv_genotype){
  
  # Unfiltered number of variants 
  stats_tab  <- data.frame(Metric = "total_unfiltered",
                           Value = nrow(pv_genotype)) 
  # Reason for filtering 
  filters_n <- table(pv_genotype$FILTER) %>%
    as.data.frame() %>%
    dplyr::rename(Metric = Var1,
                  Value = Freq) %>%
    arrange(desc(Value))
  
  stats_tab <- rbind(stats_tab,filters_n) 
  
  return(stats_tab)
}

# Function to extract genotypes for a timpoint 
get_genotype_pv <- function(PV_IDs, genotype_pv, timepoint){
  
  if (nrow(genotype_pv) > 0){
    genotype_pv <- genotype_pv %>%
      dplyr::filter(PV_ID %in% PV_IDs) 
    
    if (nrow(genotype_pv) > 0){
      genotype_pv <- genotype_pv %>%
        dplyr::select(PV_ID, AD_qual_tumour, AD_total_tumour, PVAF_qual, PVAF_total)
      colnames(genotype_pv)[-1] <- paste(timepoint, colnames(genotype_pv)[-1], sep = ".")
      
    } else {
      genotype_pv <- data.frame(PV_ID = PV_IDs,
                                AD_qual_tumour = 0, 
                                AD_total_tumour = 0, 
                                PVAF_qual = 0,
                                PVAF_total = 0)
      colnames(genotype_pv)[-1] <- paste(timepoint, colnames(genotype_pv)[-1], sep = ".")
    }
  }
  else {
    genotype_pv <- data.frame(PV_ID = PV_IDs,
                              AD_qual_tumour = 0, 
                              AD_total_tumour = 0, 
                              PVAF_qual = 0,
                              PVAF_total = 0)
    colnames(genotype_pv)[-1] <- paste(timepoint, colnames(genotype_pv)[-1], sep = ".")
  }
  
  return(genotype_pv)
}
