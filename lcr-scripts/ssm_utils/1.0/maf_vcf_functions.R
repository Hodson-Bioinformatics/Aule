##### GENERAL FUNCTIONS FOR PROCESSING MAF AND VCF FILES

### Function to load MAF
load_maf <- function(maf_files, sample_names, whitelist_gr, drivers){
  mafs <-suppressMessages(suppressWarnings(
    maf_files %>%
        map(read_tsv, col_types = cols(.default = "c"), comment = "#") %>%
        set_names(sample_names) %>%
        map(~mark_on_targets(.x, format = "maf", panel)) %>%
        map(~mark_whitelist(.x, format = "maf", whitelist_gr)) %>%
        map(~left_join(.x, drivers, by = c("Hugo_Symbol" = "Gene"))) %>%
        map(~assign_unique_ID_maf(.x)) %>%
        map(~dplyr::select(.x, any_of(c("Hugo_Symbol", "Entrez_Gene_Id",
                                        "Center", "NCBI_Build", "Chromosome", "Start_Position",
                                        "End_Position", "Strand", "Variant_Classification",
                                        "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
                                        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "t_depth",
                                        "t_ref_count", "t_alt_count", "n_depth", "n_ref_count",
                                        "n_alt_count", "on_target", "target_name", "whitelist",
                                        "Driver_Class", "vcf_pos", "VARID"))))
    ))
  return(mafs)
}

### Count variants in a maf file by type
count_by_type <- function(mafs, drivers = drivers, direct2hmrn = NULL){
  mafs %>%
    map_df(count_variants, .id = "sample") %>%
    left_join(drivers, by = c("Hugo_Symbol" = "Gene")) %>%
    mutate(group = splitvec(sample, "_", c(2), merge = "_")) %>%
    ungroup() %>%
    {if(!is.null(direct2hmrn))
        mutate(., unified_effect = direct2hmrn[Variant_Classification])
      else .}
}

### Filter variants (clinically relevant, compute rates)
filter_maf <- function(variant_types, study = "DIRECT"){

  variant_types_shared <- variant_types %>%
    dplyr::filter(!is.na(unified_effect))  %>%
    group_by(unified_effect, Hugo_Symbol, group, sample) %>%
    mutate(unified_effect_multiHit = case_when(n_distinct(unified_effect) > 1 ~ "multiHit_multiType",
                                               n_distinct(unified_effect) == 1 & sum(n) > 1 ~ "multiHit_singleType",
                                               TRUE ~ unified_effect)) %>%
    ungroup() %>%
    group_by(unified_effect_multiHit, Hugo_Symbol, group) %>%
    summarise(n_samples = n_distinct(sample),
              n_total = sum(n)) %>%
    mutate(freq = case_when(group == "Plasma" ~ n_samples/n_plasmas,
                            group == "Biopsy" ~ n_samples/n_biopsies),
           Study = study)
}

## Find on targets (format = "maf" or "vcf")
mark_on_targets <- function(variants, format = "maf", panel){
  if(format == "maf"){
    if(nrow(variants) > 0){
      variants$on_target <- FALSE
      variants$target_name <- NA
      variants_gr <- GRanges(seqnames = variants$Chromosome,
                        ranges = IRanges(start = as.numeric(variants$Start_Position),
                                         end = as.numeric(variants$End_Position)),
                        strand = "*")
      olaps <- GenomicRanges::findOverlaps(variants_gr, panel)
      variants$on_target[queryHits(olaps)] <- TRUE
      variants$target_name[queryHits(olaps)] <- panel$unified_name[subjectHits(olaps)]
    }
  } else if(format == "vcf"){
  ## TO-DO: Add VCF functionality


  }
  return(variants)
}

## Function to return stats of on target variants
on_target_stats <- function(maf){
  if(nrow(maf) > 0){
    data.frame(on_target = sum(maf$on_target), total = nrow(maf))
  }
}

## Function to count the number of variants per gene per sample in a MAF file
count_variants <- function(maf){
  if(nrow(maf) > 0){
    suppressMessages(
      maf %>%
        dplyr::filter(on_target) %>%
        group_by(Hugo_Symbol, Variant_Classification, Variant_Type) %>%
        summarise(n = n())
    )
  }
}

## Find whitelist
mark_whitelist <- function(variants, format = "maf", whitelist_gr){
  if(format == "maf"){
    if(nrow(variants) > 0){
      variants$whitelist <- FALSE
      variants_gr <- GRanges(seqnames = variants$Chromosome,
                        ranges = IRanges(start = as.numeric(variants$Start_Position),
                                         end = as.numeric(variants$End_Position)),
                        strand = "*")
      olaps <- countOverlaps(variants_gr, whitelist_gr, ignore.strand=T)
      variants$whitelist[which(olaps > 0)] <- TRUE
    }
  } else if(format == "vcf"){
  ## TO-DO: Add VCF functionality
  }
  return(variants)
}

## Metrics for individual genotypes
get_AF <- function(tidyvcf, maf){
  fix <- tidyvcf$fix
  gt <- tidyvcf$gt

  maf_to_join <- maf %>%
    mutate(REF = Reference_Allele,
           ALT = Tumor_Seq_Allele2,
           CHROM = Chromosome,
           POS = vcf_pos) %>%
    dplyr::select(target_name, Variant_Classification, Hugo_Symbol, REF, ALT, CHROM, POS, var_ID)

  integrated <- gt %>%
    left_join(fix) %>%
    dplyr::select(CHROM, POS, Indiv, gt_AF, gt_DP, REF, ALT) %>%
    separate_rows(gt_AF, ALT, sep = ",") %>%
    mutate(ALT = case_when(nchar(REF) > 1 & nchar(ALT) == 1 ~ "-",
                           nchar(REF) > 1 & nchar(ALT) > 1 ~ ALT,
                           TRUE ~ ALT),
           REF = case_when(nchar(REF) > 1 ~ substr(REF, 2, nchar(REF)),
                           TRUE ~ REF),
           gt_AF = as.numeric(gt_AF),
           Patient = splitvec(Indiv, "_", 1)) %>%
    left_join(maf_to_join, by = c("CHROM", "POS", "REF", "ALT")) %>%
    arrange(Hugo_Symbol, desc(gt_AF))

  return(integrated)
}

# Function to assign. unique variant IDs in a MAF file
assign_unique_ID_maf <- function(maf){
  maf <- maf %>%
    mutate(VARID = paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_"))
  return(maf)
}

# Find overlaping, on target variants in 2 maf files
find_overlaps_between_mafs <- function(maf1, maf2, labels){
  olaps <- data.frame(var_ID = unique(c(maf1$var_ID, maf2$var_ID))) %>%
    mutate(in_maf1 = var_ID %in% maf1$var_ID,
           in_maf2 = var_ID %in% maf2$var_ID)
  colnames(olaps) <- c("var_ID", paste0("in_", labels))
  return(olaps)
}

# Mapping between HMRN and DIRECT
rename_hmrn <- c(
  "Disruptive" = "Disruptive",
  "INFRAME INDEL" = "In_Frame_Indel",
  "Missense" = "Missense",
  "Splice" = "Splice",
  "loss of stop" = "Nonstop",
  "unknown" = "unknown"
)

direct2hmrn <- c(
  "Frame_Shift_Del" = "Disruptive",
  "Frame_Shift_Ins" = "Disruptive",
  "Nonsense_Mutation" = "Disruptive",
  "Translation_Start_Site" = "Disruptive",
  "In_Frame_Del" = "In_Frame_Indel",
  "In_Frame_Ins" = "In_Frame_Indel",
  "Missense_Mutation" = "Missense",
  "Splice_Region" = "Splice",
  "Splice_Site" = "Splice",
  "Nonstop_Mutation" = "Nonstop"
)
