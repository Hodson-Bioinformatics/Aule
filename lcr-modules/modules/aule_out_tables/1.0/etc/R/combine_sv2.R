# Author: Natasha Cutmore
# Contributor: Ilias Moutsopoulos
# Date: 08/03/2023
# Title: Join structural variant data
# Description: Join structural variant data

## If R is installed with conda, this will configure the 
## library paths to only point to the conda R package directory
local({
  conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
  rhome <- normalizePath(R.home())
  if (!is.na(conda) && startsWith(rhome, conda)) {
    .libPaths(R.home("library"), include.site = FALSE)
  }
})

# Parse snakemake options
in_csv         <- snakemake@input[["csv"]]
out_csv_all    <- snakemake@output[["csv_all"]]
out_csv_roi    <- snakemake@output[["csv_roi"]]
out_csv_simple <- snakemake@output[["csv_simple"]]
sample_info    <- snakemake@input[["sample_info"]]
log_stdout     <- snakemake@log[["stdout"]]
log_stderr     <- snakemake@log[["stderr"]]
column_name    <- snakemake@params[["join_by"]]


# Capture output
file_stderr <- file(log_stderr, open = "at")
sink(file_stderr, type = "message", append = TRUE)
sink(log_stdout, type = "output", append = TRUE)

# Load Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(bit64))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

`%!in%` <- Negate(`%in%`)

# Create, join, save table
csv_df <- tibble::tibble()
for(ii in seq_along(in_csv)){ 
  
  temp <- readr::read_csv(in_csv[ii]) %>%
    dplyr::mutate(GRIDSS.ID = NA_character_,
                  Manta.ID = NA_character_) %>%
    dplyr::mutate(sample_id = unlist(strsplit(basename(in_csv[ii]), "--"))[1], .before = 1,
                  IgCaller.Partner2_gene = as.character(IgCaller.Partner2_gene),
                  SURVIVOR.Partner1_chr = case_when(is.na(SURVIVOR.Partner1_chr) ~ NA_character_,
                                                    TRUE ~ SURVIVOR.Partner1_chr),
                  SURVIVOR.Partner1_gene = case_when(is.na(SURVIVOR.Partner1_gene) ~ NA_character_,
                                                     TRUE ~ SURVIVOR.Partner1_gene),  
                  SURVIVOR.Partner2_chr = case_when(is.na(SURVIVOR.Partner2_chr) ~ NA_character_,
                                                    TRUE ~ SURVIVOR.Partner2_chr),
                  SURVIVOR.Partner2_gene = case_when(is.na(SURVIVOR.Partner2_gene) ~ NA_character_,
                                                     TRUE ~ SURVIVOR.Partner2_gene),
                  SURVIVOR.SV_ID = case_when(is.na(SURVIVOR.SV_ID) ~ NA_character_,
                                             TRUE ~ SURVIVOR.SV_ID),
                  GRIDSS.Partner1_chr = case_when(is.na(GRIDSS.Partner1_chr) ~ NA_character_,
                                                  TRUE ~ GRIDSS.Partner1_chr),
                  GRIDSS.Partner1_gene = case_when(is.na(GRIDSS.Partner1_gene) ~ NA_character_,
                                                   TRUE ~ GRIDSS.Partner1_gene),  
                  GRIDSS.Partner2_chr = case_when(is.na(GRIDSS.Partner2_chr) ~ NA_character_,
                                                  TRUE ~ GRIDSS.Partner2_chr),
                  GRIDSS.Partner2_gene = case_when(is.na(GRIDSS.Partner2_gene) ~ NA_character_,
                                                   TRUE ~ GRIDSS.Partner2_gene),
                  GRIDSS.SV_ID = case_when(is.na(GRIDSS.SV_ID) ~ NA_character_,
                                           TRUE ~ GRIDSS.SV_ID),
                  # GRIDSS.ID = case_when(is.na(GRIDSS.ID) ~ NA_character_,
                  #                          TRUE ~ GRIDSS.ID),
                  Manta.Partner1_chr = case_when(is.na(Manta.Partner1_chr) ~ NA_character_,
                                                 TRUE ~ Manta.Partner1_chr),
                  Manta.Partner1_gene = case_when(is.na(Manta.Partner1_gene) ~ NA_character_,
                                                  TRUE ~ Manta.Partner1_gene),  
                  Manta.Partner2_chr = case_when(is.na(Manta.Partner2_chr) ~ NA_character_,
                                                 TRUE ~ Manta.Partner2_chr),
                  Manta.Partner2_gene = case_when(is.na(Manta.Partner2_gene) ~ NA_character_,
                                                  TRUE ~ Manta.Partner2_gene),
                  Manta.SV_ID = case_when(is.na(Manta.SV_ID) ~ NA_character_,
                                          TRUE ~ Manta.SV_ID),
                  IgCaller.Partner1_chr = case_when(is.na(IgCaller.Partner1_chr) ~ NA_character_,
                                                    TRUE ~ IgCaller.Partner1_chr),
                  IgCaller.Partner1_gene = case_when(is.na(IgCaller.Partner1_gene) ~ NA_character_,
                                                     TRUE ~ IgCaller.Partner1_gene),  
                  IgCaller.Partner2_chr = case_when(is.na(IgCaller.Partner2_chr) ~ NA_character_,
                                                    TRUE ~ IgCaller.Partner2_chr),
                  IgCaller.Partner2_gene = case_when(is.na(IgCaller.Partner2_gene) ~ NA_character_,
                                                     TRUE ~ IgCaller.Partner2_gene),
                  IgCaller.SV_ID = case_when(is.na(IgCaller.SV_ID) ~ NA_character_,
                                             TRUE ~ IgCaller.SV_ID),
                  SV_ID = case_when(is.na(SV_ID) ~ NA_character_,
                                    TRUE ~ SV_ID)) %>%
    dplyr::filter(!is.na(SV_ID) | !is.na(IgCaller.SV_ID)) %>%
    distinct()
  
  if (nrow(csv_df) == 0) {
    csv_df <- temp
  } else {
    csv_df <- csv_df %>%
      full_join(temp)  
  }
}


#Create output tables
clin_sig_chr <- c("chr2", "chr3", "chr8", "chr14", "chr18", "chr22")

fill_cols <- c("GRIDSS.Partner1_chr", "GRIDSS.Partner1_gene",
               "GRIDSS.Partner2_chr", "GRIDSS.Partner2_gene",
               "Manta.Partner1_chr", "Manta.Partner1_gene",
               "Manta.Partner2_chr", "Manta.Partner2_gene",
               "IgCaller.Partner1_chr", "IgCaller.Partner1_gene",
               "IgCaller.Partner2_chr", "IgCaller.Partner2_gene")

sample_df <- readr::read_csv(sample_info) %>%
  dplyr::filter(direct_timepoint == "B") %>%
  dplyr::select(sample_id, patient_id, direct_type, direct_timepoint)

by_all.df <- csv_df %>%
  dplyr::left_join(sample_df, by = "sample_id") %>%
  dplyr::relocate(patient_id, .after = sample_id) %>%
  dplyr::relocate(direct_type, .after = patient_id) %>%
  dplyr::select(-c(ends_with("pos"),
                   starts_with("SURVIVOR"),
                   starts_with("Found"),
                   ends_with("SV_ID"),
                   ends_with("VAF"),
                   "GRIDSS.ID", "Manta.ID")) %>%
  distinct() %>%
  dplyr::mutate(IgCaller.Partner1_gene = case_when(grepl("^IG", IgCaller.Partner1_gene)
                                                   ~ substr(IgCaller.Partner1_gene, 1, 3),
                                                   TRUE ~ IgCaller.Partner1_gene),
                IgCaller.Partner1_gene = case_when(IgCaller.Partner1_gene == "LOC105372161"
                                                   ~ "BCL2", TRUE ~ IgCaller.Partner1_gene),
                IgCaller.Partner2_gene = case_when(grepl("^IG", IgCaller.Partner2_gene)
                                                   ~ substr(IgCaller.Partner2_gene, 1, 3),
                                                   TRUE ~ IgCaller.Partner2_gene),
                IgCaller.Partner2_gene = case_when(IgCaller.Partner2_gene == "LOC105372161"
                                                   ~ "BCL2", TRUE ~ IgCaller.Partner2_gene),
                IgCaller.Partner1_gene = case_when(IgCaller.Partner1_chr %!in% clin_sig_chr
                                                   ~ NA_character_, TRUE ~ IgCaller.Partner1_gene),
                IgCaller.Partner2_gene = case_when(IgCaller.Partner1_chr %!in% clin_sig_chr
                                                   ~ NA_character_, TRUE ~ IgCaller.Partner2_gene),
                IgCaller.Partner2_chr = case_when(IgCaller.Partner1_chr %!in% clin_sig_chr
                                                  ~ NA_character_, TRUE ~ IgCaller.Partner2_chr),
                IgCaller.Partner1_chr = case_when(IgCaller.Partner1_chr %!in% clin_sig_chr
                                                  ~ NA_character_, TRUE ~ IgCaller.Partner1_chr),
                IgCaller.Partner1_gene = case_when(IgCaller.Partner2_chr %!in% clin_sig_chr
                                                   ~ NA_character_, TRUE ~ IgCaller.Partner1_gene),
                IgCaller.Partner2_gene = case_when(IgCaller.Partner2_chr %!in% clin_sig_chr
                                                   ~ NA_character_, TRUE ~ IgCaller.Partner2_gene),
                IgCaller.Partner1_chr = case_when(IgCaller.Partner2_chr %!in% clin_sig_chr
                                                  ~ NA_character_, TRUE ~ IgCaller.Partner1_chr),
                IgCaller.Partner2_chr = case_when(IgCaller.Partner2_chr %!in% clin_sig_chr
                                                  ~ NA_character_, TRUE ~ IgCaller.Partner2_chr)) %>% 
  dplyr::group_by(sample_id, patient_id, direct_type) %>%
  fill(all_of(fill_cols), .direction = "downup") %>%
  ungroup() %>%
  distinct() %>%
  dplyr::relocate(patient_id, .after = sample_id) %>%
  dplyr::relocate(direct_type, .after = sample_id) 
write.csv(by_all.df, out_csv_all, row.names = FALSE)



by_sample.df <- expand.grid(
  patient_id = sample_df$patient_id,
  direct_type = by_all.df$direct_type
) %>%
  dplyr::left_join(sample_df) %>%
  dplyr::left_join(by_all.df) %>%
  distinct() %>%
  dplyr::mutate(BCL2 = case_when(
                                  (GRIDSS.Partner1_gene == "BCL2" |
                                     GRIDSS.Partner2_gene == "BCL2" |
                                     Manta.Partner1_gene == "BCL2" |
                                     Manta.Partner2_gene == "BCL2" |
                                     IgCaller.Partner1_gene == "BCL2" |
                                     IgCaller.Partner2_gene == "BCL2") ~1, TRUE ~ 0),
                BCL6 = case_when(
                                  (GRIDSS.Partner1_gene == "BCL6" |
                                     GRIDSS.Partner2_gene == "BCL6" |
                                     Manta.Partner1_gene == "BCL6" |
                                     Manta.Partner2_gene == "BCL6" |
                                     IgCaller.Partner1_gene == "BCL6" |
                                     IgCaller.Partner2_gene == "BCL6") ~1, TRUE ~ 0),
                MYC = case_when(
                                  (GRIDSS.Partner1_gene == "MYC" |
                                     GRIDSS.Partner2_gene == "MYC" |
                                     Manta.Partner1_gene == "MYC" |
                                     Manta.Partner2_gene == "MYC" |
                                     IgCaller.Partner1_gene == "MYC" |
                                     IgCaller.Partner2_gene == "MYC") ~1, TRUE ~ 0)) %>%
  dplyr::mutate_if(is.numeric, ~replace_na(., 0)) %>%
  distinct() %>%
  dplyr::relocate(sample_id, .before = 1) %>%
  dplyr::select(sample_id, patient_id, direct_type, BCL2, BCL6, MYC) %>%
write.csv(out_csv_roi, row.names = FALSE)

by_patient.df <- expand.grid(
  patient_id = sample_df$patient_id,
  direct_type = by_all.df$direct_type
) %>%
  dplyr::left_join(by_all.df) %>%
  distinct() %>%
  dplyr::mutate(BCL2 = case_when(
                                  (GRIDSS.Partner1_gene == "BCL2" |
                                   GRIDSS.Partner2_gene == "BCL2" |
                                   Manta.Partner1_gene == "BCL2" |
                                   Manta.Partner2_gene == "BCL2" |
                                   IgCaller.Partner1_gene == "BCL2" |
                                   IgCaller.Partner2_gene == "BCL2") ~1, TRUE ~ 0),
                BCL6 = case_when(
                                  (GRIDSS.Partner1_gene == "BCL6" |
                                   GRIDSS.Partner2_gene == "BCL6" |
                                   Manta.Partner1_gene == "BCL6" |
                                   Manta.Partner2_gene == "BCL6" |
                                   IgCaller.Partner1_gene == "BCL6" |
                                   IgCaller.Partner2_gene == "BCL6") ~1, TRUE ~ 0),
                MYC = case_when(
                                  (GRIDSS.Partner1_gene == "MYC" |
                                   GRIDSS.Partner2_gene == "MYC" |
                                   Manta.Partner1_gene == "MYC" |
                                   Manta.Partner2_gene == "MYC" |
                                   IgCaller.Partner1_gene == "MYC" |
                                   IgCaller.Partner2_gene == "MYC") ~1, TRUE ~ 0)) %>%
  dplyr::select(patient_id, BCL2, BCL6, MYC) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::summarise(BCL2 = sum(BCL2),
                   BCL6 = sum(BCL6),
                   MYC = sum(MYC)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(across(where(is.numeric), .fns = ~replace(., . >= 1 , 1))) %>%
  distinct() %>%
write.csv(out_csv_simple, row.names = FALSE)

# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)

save.image()
