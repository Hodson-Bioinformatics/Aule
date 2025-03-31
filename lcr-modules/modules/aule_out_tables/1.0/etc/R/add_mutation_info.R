# Author: Ilias Moutsopoulos
# Contributor: Natasha Cutmore
# Date: 02/03/2023
# Title: Handle mutation data
# Description: Concatenate mutation data and add whitelist/driver information

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
vcfs           <- snakemake@input[["vcfs"]]
mafs           <- snakemake@input[["mafs"]]
whitelist_path <- snakemake@input[["whitelist"]]
drivers_path   <- snakemake@input[["drivers"]]
out_csv        <- snakemake@output[["csv"]]
out_csv_clean  <- snakemake@output[["csv_clean"]]
vaf_filter     <- snakemake@params[["vaf_filter"]]
shm_distance   <- snakemake@params[["shm_distance"]]
log_stdout     <- snakemake@log[["stdout"]]
log_stderr     <- snakemake@log[["stderr"]]

# Capture output
file_stderr <- file(log_stderr, open = "at")
sink(file_stderr, type = "message", append = TRUE)
sink(log_stdout, type = "output", append = TRUE)

# Load Libraries
library(readr)
library(tibble)
library(dplyr)
library(purrr)
library(vcfR)
library(GenomicRanges)
library(biomaRt)

`%!in%` <- Negate(`%in%`)

# Read whitelist and drivers
whitelist <- readr::read_csv(whitelist_path)
drivers_list <- readr::read_delim(drivers_path)

# Define helpers
TSG_genes <- drivers_list$Gene[drivers_list$Driver_Class == "TSG"]
TSG_consequences <- c("frameshift_variant", "stop_lost", "stop_gained",
                      "start_lost", "splice_donor_variant",
                      "splice_acceptor_variant") # , "start_gained", "protein_altering_variant"
SHM_genes <- c(drivers_list$Gene[drivers_list$Driver_Class == "SHM"], "BCL6" )
oncogenes <- drivers_list$Gene[drivers_list$Driver_Class %in% c("Oncogene", "Oncogenes")]
noncoding_consequences <- c("5'UTR", "5'Flank", "3'UTR", "3'Flank", "Intron")
between_coords <- purrr::pmap(
  .l = list(
    drivers_list$Chr,
    drivers_list$transcription_start_site - shm_distance,
    drivers_list$transcription_start_site + shm_distance
  ),
  .f = ~ paste0(..1, "_", ..2:..3)
) %>%
  unlist() %>%
  unique()

# Read, tidy and combine all mafs
maf_df <- tibble()
for(ii in seq_along(mafs)){
  temp <- readr::read_delim(mafs[ii]) %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position,
                  Reference_Allele, Tumor_Seq_Allele2, Exon_Number,
                  HGVSc, HGVSp, Variant_Classification, Variant_Type, Consequence, # added HGVSc for Seshat joining
                  Entrez_Gene_Id, Center, NCBI_Build, Transcript_ID) %>%
    dplyr::mutate(sample_id = unlist(strsplit(basename(mafs[ii]), "--"))[1], .before = 1) %>%
    dplyr::mutate(
      HG38_uniqueID = paste(Chromosome, Start_Position,
                            Reference_Allele, Tumor_Seq_Allele2,
                            sep = "_"),
      Sample_mut_ID = paste0(sample_id, ":", HG38_uniqueID),
      Start_Position = as.numeric(Start_Position),
      Tumor_Seq_Allele2 = as.character(Tumor_Seq_Allele2),
      Tumor_Seq_Allele2 = ifelse(Tumor_Seq_Allele2 == "TRUE", "T", Tumor_Seq_Allele2), # fix TRUE quirk
      Entrez_Gene_Id =  NULL
    )

  maf_df <- rbind(maf_df, temp)
}

# Get ENTREZ ID as vcf2maf always outputs 0!!!
biomart_hg38 <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "entrezgene_id",
                 "transcription_start_site", "strand"),
  filters = "ensembl_transcript_id",
  values = maf_df$Transcript_ID,
  mart = biomaRt::useEnsembl(
    biomart="ensembl",
    dataset="hsapiens_gene_ensembl",
    host = "https://jul2019.archive.ensembl.org"
  )
) %>%
  dplyr::distinct() %>%
  dplyr::rename(Entrez_Gene_Id = entrezgene_id, TSS = transcription_start_site)

# Read, tidy and combine all vcfs
vcf_df <- tibble()
for(ii in seq_along(vcfs)){
  vcf_in <- vcfR::read.vcfR(vcfs[ii])
  if(nrow(vcf_in) == 0){
    temp <- tibble::tibble()
  }else{
    vcf_in_tidy <- vcfR::vcfR2tidy(vcf_in)

    temp <- vcf_in_tidy$gt %>%
      dplyr::filter(Indiv == "TUMOR") %>%
      dplyr::select(-c("POS")) %>%
      cbind(vcf_in_tidy$fix) %>%
      dplyr::select(-c(ChromKey)) %>%
      dplyr::filter(FILTER == "PASS") %>%
      dplyr::select(CHROM, POS, REF, ALT, gt_AF) %>%
      dplyr::filter(gt_AF > vaf_filter) %>%
      dplyr::mutate(sample_id = unlist(strsplit(basename(vcfs[ii]), "--"))[1], .before = 1) %>%
      dplyr::mutate(HG38_uniqueID = paste(CHROM, POS, REF, ALT, sep = "_"),
                    Sample_mut_ID = paste0(sample_id, ":", HG38_uniqueID)) %>%
      tibble::as_tibble()
  }

  vcf_df <- rbind(vcf_df, temp)
}

# Generate output df
df <- vcf_df %>%
  # Align deletions and insertions encoding between vcf and maf
  dplyr::mutate( # deletions
    vcf.start = ifelse(Sample_mut_ID %!in% maf_df$Sample_mut_ID & nchar(ALT) < nchar(REF), POS + 1, POS),
    vcf.ALT = ifelse(Sample_mut_ID %!in% maf_df$Sample_mut_ID & nchar(ALT) < nchar(REF), "-", ALT),
    vcf.REF = ifelse(Sample_mut_ID %!in% maf_df$Sample_mut_ID & nchar(ALT) < nchar(REF), stringr::str_sub(REF, 2), REF),
    HG38_uniqueID = ifelse(vcf.start != POS, paste0(CHROM, "_", vcf.start, "_", vcf.REF, "_", vcf.ALT), HG38_uniqueID),
    Sample_mut_ID = ifelse(vcf.start != POS, paste0(sample_id, ":", CHROM, "_", vcf.start, "_", vcf.REF, "_", vcf.ALT), Sample_mut_ID)
  ) %>%
  dplyr::mutate( # insertions
    vcf.REF = ifelse(Sample_mut_ID %!in% maf_df$Sample_mut_ID & nchar(ALT) > nchar(REF), "-", REF),
    vcf.ALT = ifelse(Sample_mut_ID %!in% maf_df$Sample_mut_ID & nchar(ALT) > nchar(REF), stringr::str_sub(ALT, 2), ALT),
    HG38_unqiueID = ifelse(vcf.ALT != ALT, paste0(CHROM, "_", vcf.start, "_", vcf.REF, "_", vcf.ALT), HG38_uniqueID),
    Sample_mut_ID = ifelse(vcf.ALT != ALT, paste0(sample_id, ":", CHROM, "_", vcf.start, "_", vcf.REF, "_", vcf.ALT), Sample_mut_ID)
  ) %>%
  # Join all information
  dplyr::left_join(maf_df, by = c("Sample_mut_ID", "sample_id", "HG38_uniqueID")) %>%
  dplyr::left_join(biomart_hg38, by = c("Transcript_ID" = "ensembl_transcript_id")) %>%
  # Add rules for classifying mutations
  dplyr::mutate(
    Whitelist = ifelse(HG38_uniqueID %in% whitelist$HG38_uniqueID, 1, 0),
    Driver_Mutation = case_when(
      (Hugo_Symbol %in% TSG_genes & Consequence %in% TSG_consequences) ~ 1,
      (Hugo_Symbol %in% TSG_genes & grepl("^missense", Consequence) & Whitelist == 1) ~ 1,
      (Hugo_Symbol %in% TSG_genes & grepl("^In_Frame", Variant_Classification) & Whitelist == 1) ~ 1,
      (Hugo_Symbol %in% TSG_genes & grepl("Splice_Site", Variant_Classification)) ~ 1,
      (Hugo_Symbol %in% oncogenes & grepl("^missense", Consequence) & Whitelist == 1) ~ 1,
      (Hugo_Symbol %in% oncogenes & grepl("^Splice", Variant_Classification) & Whitelist == 1) ~ 1,
      (Hugo_Symbol %in% oncogenes & grepl("^In_Frame", Variant_Classification) & Whitelist == 1) ~ 1,
      (Hugo_Symbol == "NOTCH2" & Exon_Number == "34/34" & Consequence %in% TSG_consequences) ~ 1, # CHROM == "chr1" & dplyr::between(POS, 119911553, 119916694)
      (Hugo_Symbol == "NOTCH1" & Exon_Number == "34/34" & Consequence %in% TSG_consequences) ~ 1, # CHROM == "chr9" & dplyr::between(POS, 136494433, 136497558)
      (Hugo_Symbol %in% c("NOTCH1", "NFKBIZ") & Variant_Classification == "3'UTR") ~ 1,
      (Hugo_Symbol == "TP53" & (Consequence %in% TSG_consequences | Whitelist == 1)) ~ 1,
      TRUE ~ 0
    ),
    SHM_TSS_proximal_Mutation =  ifelse((Hugo_Symbol %in% SHM_genes & paste0(CHROM, "_", POS) %in% between_coords), 1, 0), # Added BCL6 as this is not listed as SHM gene in driver list
    SHM_Noncoding_Mutation = ifelse((Hugo_Symbol %in% SHM_genes & Variant_Classification %in% noncoding_consequences), 1, 0),
    SHM_Coding_Mutation = ifelse((Hugo_Symbol %in% SHM_genes & Driver_Mutation == 1), 1, 0),
    Hotspot_Mutation = case_when(
      (Hugo_Symbol == "MYD88" & grepl("p.Leu273Pro", HGVSp)) ~ 1,
      (Hugo_Symbol == "EZH2" & grepl("^p.Tyr646", HGVSp)) ~ 1,
      (Hugo_Symbol == "POU2F2" & grepl("^p.Thr239", HGVSp)) ~ 1,
      (Hugo_Symbol == "MYC" & grepl("^p.Thr58", HGVSp)) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  # Remove duplicate  columns
  dplyr::select(-c(dplyr::starts_with("vcf"), Chromosome, Start_Position,
                   Reference_Allele, Tumor_Seq_Allele2, Sample_mut_ID))

# Save in format friendly for maftools visualisation
write_csv(df, out_csv)

# Save in human friendly version
df %>%
  dplyr::rename(VAF = gt_AF, Gene = Hugo_Symbol, AminoAcid_change = HGVSp) %>%
  dplyr::select(-c(HG38_uniqueID, Variant_Type, #Entrez_Gene_Id, ## needed in LG!
                   Center, NCBI_Build, Transcript_ID)) %>%
  write_csv(out_csv_clean)

# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)
