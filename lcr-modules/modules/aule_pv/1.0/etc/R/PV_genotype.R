##### Script ot perform basic filtering (R1/R2 concordance, MQ, BQ)
suppressPackageStartupMessages(library(tidyverse))
snakemake@source("PV_utilis.R")

genotype_pv <- function(pv_annotated, pv_type, max_matched_germline_support, max_cohort_germline_support, max_gnomad_support, max_pon_snv_support, max_pon_pv_support, min_AF, min_DP){

  # Load 
  pv_tab_raw <- read_delim(pv_annotated) 
  
  # Filter sample based on provided criteria 
  if (pv_type == "doublets"){
    pv_tab <- pv_tab_raw %>%
      ungroup() %>%
      mutate(
        total_DP = pmin(DP_pv1, DP_pv2),
        total_DP = pmin(DP_pv1, DP_pv2),
        PVAF_qual = AD_qual_tumour/total_DP,
        PVAF_total = AD_total_tumour/total_DP,
        germline_total_DP = pmin(germline_DP_pv1, germline_DP_pv2),
        PON_PV_Nsamples = as.numeric(PON_PV_Nsamples),
        # Germline support
        FILTER = case_when(germline_AD_pv1 > max_matched_germline_support | 
                           germline_AD_pv2 > max_matched_germline_support | 
                           AD_total_germline > max_matched_germline_support ~ "germline",
                           TRUE ~ ""),
        # Max allele frequency in gnomad 
        FILTER = case_when(gnomad_pv1 > max_gnomad_support | gnomad_pv2 > max_gnomad_support ~ paste(FILTER, "gnomad", sep = ";"),
                           TRUE ~ FILTER),
        # Max PON SNV support
        FILTER = case_when(pon_pv1 > max_pon_snv_support | pon_pv2 > max_pon_snv_support ~ paste(FILTER, "pon_snv", sep = ";"),
                           TRUE ~ FILTER),
        # TBA:Max PON PV support
        FILTER = case_when(PON_PV_Nsamples > max_pon_pv_support ~ paste(FILTER, "pon_pv", sep = ";"),
                           TRUE ~ FILTER),
        # min AF
        FILTER = case_when(PVAF_qual < min_AF ~ paste(FILTER, "low_pvaf", sep = ";"),
                           TRUE ~ FILTER),
        # min DP
        FILTER = case_when(total_DP < min_DP ~ paste(FILTER, "low_dp", sep = ";"),
                           TRUE ~ FILTER),
        FILTER = case_when(FILTER == "" ~ "PASS",
                           grepl("^;", FILTER) ~ substr(FILTER, 2, nchar(FILTER)),
                           TRUE ~ FILTER))
      # TBA: Max cohort germline support 
      #dplyr::filter(germline_cohort_pv1 <= max_cohort_germline_support & germline_cohort_pv2 <= max_cohort_germline_support) %>% 
  }
  else if (pv_type == "triplets"){
    pv_tab <- pv_tab_raw %>%
      ungroup() %>%
      mutate(
        total_DP = pmin(DP_pv1, DP_pv2, DP_pv3),
        total_DP = pmin(DP_pv1, DP_pv2, DP_pv3),
        PVAF_qual = AD_qual_tumour/total_DP,
        PVAF_total = AD_total_tumour/total_DP,
        germline_total_DP = pmin(germline_DP_pv1, germline_DP_pv2, germline_DP_pv3),
        # Germline support
        FILTER = case_when(germline_AD_pv1 > max_matched_germline_support | 
                           germline_AD_pv2 > max_matched_germline_support | 
                             germline_AD_pv3 > max_matched_germline_support | 
                             AD_total_germline > max_matched_germline_support ~ "germline",
                           TRUE ~ ""),
        # Max allele frequency in gnomad 
        FILTER = case_when(gnomad_pv1 > max_gnomad_support | gnomad_pv2 > max_gnomad_support | gnomad_pv3 > max_gnomad_support ~ paste(FILTER, "gnomad", sep = ";"),
                           TRUE ~ FILTER),
        # Max PON SNV support
        FILTER = case_when(pon_pv1 > max_pon_snv_support | pon_pv2 > max_pon_snv_support  | pon_pv3 > max_pon_snv_support  ~ paste(FILTER, "pon_snv", sep = ";"),
                           TRUE ~ FILTER),
        # TBA:Max PON PV support
        FILTER = case_when(PON_PV_Nsamples > max_pon_pv_support ~ paste(FILTER, "pon_pv", sep = ";"),
                           TRUE ~ FILTER),
        # min AF
        FILTER = case_when(PVAF_qual < min_AF ~ paste(FILTER, "low_pvaf", sep = ";"),
                           TRUE ~ FILTER),
        # min DP
        FILTER = case_when(total_DP < min_DP ~ paste(FILTER, "low_dp", sep = ";"),
                           TRUE ~ FILTER),
        FILTER = case_when(FILTER == "" ~ "PASS",
                           grepl("^;", FILTER) ~ substr(FILTER, 2, nchar(FILTER)),
                           TRUE ~ FILTER))
  }
 
  return(list(pv_tab = pv_tab))
}

## Load params

pv_annotated  <- snakemake@input[["pv_annotated"]]
pv_genotyped  <- snakemake@output[["pv_genotyped"]]
pv_type       <- snakemake@params[["pv_type"]]
max_matched_germline_support <- snakemake@params[["max_matched_germline_support"]]
max_cohort_germline_support  <- snakemake@params[["max_cohort_germline_support"]]
max_gnomad_support           <- snakemake@params[["max_gnomad_support"]]
max_pon_snv_support          <- snakemake@params[["max_pon_snv_support"]]
max_pon_pv_support           <- snakemake@params[["max_pon_pv_support"]]
min_AF                       <- snakemake@params[["min_AF"]]
min_DP                       <- snakemake@params[["min_DP"]]


# Collapse PVs
genotyped <- genotype_pv(pv_annotated = pv_annotated, 
                         pv_type = pv_type,
                         max_matched_germline_support = max_matched_germline_support, 
                         max_cohort_germline_support = max_cohort_germline_support, 
                         max_gnomad_support = max_gnomad_support,
                         max_pon_snv_support = max_pon_snv_support,
                         max_pon_pv_support = max_pon_pv_support,
                         min_AF = min_AF,
                         min_DP = min_DP)

write_delim(genotyped$pv_tab, file = pv_genotyped, delim = "\t")
