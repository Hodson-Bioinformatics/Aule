# Author: Ilias Moutsopoulos
# Contributor: Natasha Cutmore, Joanna Krupka
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
mafs           <- snakemake@input[["mafs"]]
whitelist_path <- snakemake@input[["whitelist"]]
drivers_path   <- snakemake@input[["drivers"]]
panel_file     <- snakemake@input[["panel"]]
blacklisted_coord_file <- snakemake@input[["blacklisted_coord"]]
out_csv        <- snakemake@output[["csv"]]
filtering      <- snakemake@params[["filtering"]]
shm_distance   <- snakemake@params[["shm_distance"]]
min_pon        <- snakemake@params[["min_pon"]]
min_germline   <- snakemake@params[["min_germline"]]
max_error_prob <- snakemake@params[["max_error_prob"]]
min_umi        <- snakemake@params[["min_umi"]]
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

# Function to compute the likelihood of error
compute_error_likelihood <- function(n, k, p) {
    # n: total depth of sequencing
    # k: minimum number of errors
    # p: probability of sequencing error

    # Compute the cumulative probability for X < k
    cumulative_prob <- pbinom(k - 1, n, p)

    # Probability of observing at least k errors
    at_least_k_prob <- 1 - cumulative_prob

    return(at_least_k_prob)
}

# Function to test the significance of strand balance between ref and alt supporting reads
test_balanced_strand_support <- function(ref_fwd, ref_rev, mut_fwd, mut_rev) {

    # Create a 2x2 contingency table
    contingency_table <- matrix(
        c(ref_fwd, ref_rev, mut_fwd, mut_rev),
        nrow = 2,
        byrow = T
    )

    # Perform Fisher's Exact Test
    test_result <- fisher.test(contingency_table)

    # Return the p-value
    return(test_result$p.value)
}

# Read whitelist and drivers
message("Loading metadata")
whitelist <- readr::read_delim(whitelist_path, show_col_types = FALSE)
drivers_list <- readr::read_delim(drivers_path, show_col_types = FALSE)
panel <- read_delim(panel_file)
panel_gr <- GRanges(panel) + 2
blacklisted_coord <- read_csv(blacklisted_coord_file)
blacklisted_coord_gr <- GRanges(seqnames = blacklisted_coord$chr,
                                IRanges(start =  blacklisted_coord$pos,
                                        end = blacklisted_coord$pos),
                                strand = "*")

# Define helpers
TSG_genes <- drivers_list$Gene[drivers_list$Driver_Class == "TSG"]
TSG_consequences <- c("frameshift_variant", "stop_lost", "stop_gained",
                      "start_lost", "splice_donor_variant",
                      "splice_acceptor_variant")
SHM_genes <- drivers_list$Gene[drivers_list$Driver_Class == "SHM"]
oncogenes <- drivers_list$Gene[drivers_list$Driver_Class %in% c("Oncogene", "Oncogenes")]
noncoding_consequences <- c("5'UTR", "5'Flank", "3'UTR", "3'Flank", "Intron", "RNA", "Silent", "IGR")
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
        dplyr::mutate(sample_id = unlist(strsplit(basename(mafs[ii]), "--"))[1], .before = 1,
                      POSID = paste(Chromosome, Start_Position, sep = "_"),
                      VARID = paste(Chromosome, Start_Position,
                                    Reference_Allele, Tumor_Seq_Allele2,
                                    sep = "_"),
                      Sample_VARID = paste0(sample_id, ":", VARID),
                      Start_Position = as.numeric(Start_Position),
                      Tumor_Seq_Allele2 = as.character(Tumor_Seq_Allele2),
                      Tumor_Seq_Allele2 = ifelse(Tumor_Seq_Allele2 == "TRUE", "T", Tumor_Seq_Allele2), # fix TRUE quirk
                      Entrez_Gene_Id =  NULL
        ) %>%
        # Remove empty columns
        dplyr::select(-dbSNP_Val_Status, -Tumor_Validation_Allele1, -Tumor_Validation_Allele2,
                      -Match_Norm_Validation_Allele1, -Match_Norm_Validation_Allele2,
                      -Verification_Status, -Validation_Status, -Mutation_Status,
                      -Sequencing_Phase, -Sequence_Source, -Validation_Method,
                      -Score, -BAM_File, -Sequencer, -Tumor_Sample_UUID, -Matched_Norm_Sample_UUID,
                      -PUBMED, -MOTIF_NAME, -MOTIF_POS, -HIGH_INF_POS, -MOTIF_SCORE_CHANGE,
                      -HGVS_OFFSET, -MINIMISED, -vcf_id)

    maf_df <- rbind(maf_df, temp)
}

# Get ENTREZ ID as vcf2maf always outputs 0!!!
message("Loading Biomart annotations")
biomart_hg38 <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "entrezgene_id",
                   "transcription_start_site", "strand"),
    filters = "ensembl_transcript_id",
    values = unique(maf_df$Transcript_ID),
    mart = biomaRt::useEnsembl(
        biomart="ensembl",
        dataset="hsapiens_gene_ensembl",
        version=105
    ))  %>%
    dplyr::distinct() %>%
    group_by(ensembl_transcript_id) %>%
    dplyr::summarise(entrezgene_id = min(entrezgene_id),
              transcription_start_site = transcription_start_site[which(entrezgene_id == min(entrezgene_id))],
              strand = unique(strand)[which(entrezgene_id == min(entrezgene_id))]) %>%
    dplyr::rename(Entrez_Gene_Id = entrezgene_id,
                  TSS = transcription_start_site)



# Identify off target regions
loci <- GRanges(seqnames = maf_df$Chromosome,
                    ranges = IRanges(start = maf_df$Start_Position,
                                     end   = maf_df$End_Position),
                    strand = "*")
olaps <- countOverlaps(loci, panel_gr)
maf_df$on_target <- FALSE
maf_df$on_target[olaps > 0] <- TRUE

print(filtering)
if (filtering == "True"){
    # Identify blacklisted loci
    olaps <- countOverlaps(loci, blacklisted_coord_gr)
    maf_df$blacklisted <- FALSE
    maf_df$blacklisted[olaps > 0] <- TRUE

    # Filter MAF (Plasma only samples)
    maf_df <- maf_df %>%
        mutate(error_prob = compute_error_likelihood(t_depth, t_alt_count, directional_error_rate),
               variant_in_3plus = ifelse(is.na(variant_in_3plus), 0, variant_in_3plus)) %>%
        rowwise() %>%
        # Pre-filetring computations
        mutate(umi_family_reference_nontrivial_revstrand = umi_family_reference_nontrivial - umi_family_reference_nontrivial_fwdstrand,
               umi_family_supporting_nontrivial_revstrand = umi_family_supporting_nontrivial - umi_family_supporting_nontrivial_fwdstrand,
               strand_balance = test_balanced_strand_support(umi_family_reference_nontrivial_fwdstrand,
                                                             umi_family_reference_nontrivial_revstrand,
                                                             umi_family_supporting_nontrivial_fwdstrand,
                                                             umi_family_supporting_nontrivial_revstrand)) %>%
        ungroup() %>%
        # Apply filters
        mutate(FILTER = ifelse(variant_in_3plus == 1, "PASS", ""),
               # PON, regardless of SLMS3 status
               FILTER = case_when(FILTER == "PASS" & PON_AF > min_pon ~ "pon;",
                                  FILTER != "PASS" & PON_AF > min_pon ~ paste0(FILTER, "pon", sep = ";"),
                                  TRUE ~ FILTER),
               # Off target, regardless of SLMS3 status
               FILTER = case_when(FILTER == "PASS" & !on_target ~ "off_target;",
                                  FILTER != "PASS" & !on_target ~ paste0(FILTER, "off_target", sep = ";"),
                                  TRUE ~ FILTER),
               # # Blacklisting, regardless of SLMS3 status
               FILTER = case_when(FILTER == "PASS" & blacklisted ~ "blacklisted;",
                                  FILTER != "PASS" & blacklisted ~ paste0(FILTER, "blacklisted", sep = ";"),
                                  TRUE ~ FILTER),
               # Max germline support for SLMS2 variants
               FILTER = case_when(FILTER != "PASS" & n_alt_count > min_germline ~ paste0(FILTER, "germline", sep = ";"),
                                  TRUE ~ FILTER),
               # Probability of error
               FILTER = case_when(FILTER != "PASS" & error_prob > max_error_prob ~ paste0(FILTER, "errorProb", sep = ";"),
                                  TRUE ~ FILTER),
               # UMI support
               FILTER = case_when(FILTER != "PASS" & umi_family_supporting_nontrivial < min_umi ~ paste0(FILTER, "umi", sep = ";"),
                                  TRUE ~ FILTER),
               #Balanced strand support
               FILTER = case_when(FILTER != "PASS" & strand_balance < max_sb ~ paste0(FILTER, "strandBias", sep = ";"),
                                  TRUE ~ FILTER),
               FILTER = ifelse(FILTER == "", "PASS", FILTER)
        ) %>%
        dplyr::relocate(VAF, FILTER, .after = "Tumor_Seq_Allele2")
}

# Generate output df
message("Annotating MAF with Driver")
df <- maf_df %>%
    ungroup() %>%
        # Join all information
        dplyr::left_join(biomart_hg38, by = c("Transcript_ID" = "ensembl_transcript_id")) %>%
        # Add rules for classifying mutations
        dplyr::mutate(
            Whitelist = ifelse(VARID %in% whitelist$HG38_uniqueID, 1, 0),
            Driver_Mutation_v1 = case_when(
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
            Driver_Mutation_v2 = case_when(
                (Hugo_Symbol %in% TSG_genes & Consequence %in% TSG_consequences) ~ 1,
                (Hugo_Symbol %in% TSG_genes & grepl("Splice_Site", Variant_Classification)) ~ 1,
                (Whitelist == 1) ~ 1,
                (Hugo_Symbol == "NOTCH2" & Exon_Number == "34/34" & Consequence %in% TSG_consequences) ~ 1, # CHROM == "chr1" & dplyr::between(POS, 119911553, 119916694)
                (Hugo_Symbol == "NOTCH1" & Exon_Number == "34/34" & Consequence %in% TSG_consequences) ~ 1, # CHROM == "chr9" & dplyr::between(POS, 136494433, 136497558)
                (Hugo_Symbol %in% c("NOTCH1", "NFKBIZ") & Variant_Classification == "3'UTR") ~ 1,
                TRUE ~ 0
            ),
            SHM_TSS_proximal_Mutation =  ifelse((Hugo_Symbol %in% SHM_genes & POSID %in% between_coords), 1, 0), # Added BCL6 as this is not listed as SHM gene in driver list
            SHM_Noncoding_Mutation = ifelse((Hugo_Symbol %in% SHM_genes & Variant_Classification %in% noncoding_consequences), 1, 0),
            SHM_Coding_Mutation = ifelse((Hugo_Symbol %in% SHM_genes & Driver_Mutation_v1 == 1), 1, 0),
            Hotspot_Mutation = case_when(
                (Hugo_Symbol == "MYD88" & grepl("p.Leu252Pro", HGVSp)) ~ 1,
                (Hugo_Symbol == "EZH2" & grepl("^p.Tyr646", HGVSp)) ~ 1,
                (Hugo_Symbol == "POU2F2" & grepl("^p.Thr261", HGVSp)) ~ 1,
                (Hugo_Symbol == "MYC" & grepl("^p.Thr58", HGVSp)) ~ 1,
                TRUE ~ 0
            )
        ) %>%
        distinct()

# Check duplicated rows
message("Checking for duplicates")
duplicated <- table(df$Sample_VARID) %>%
    as.data.frame() %>%
    dplyr::filter(Freq > 1)
df_dup <- df %>%
    dplyr::filter(Sample_VARID %in% duplicated$Var1)
if(nrow(duplicated) > 0){
    stop("Duplicated variants in the final output found, check the add_mutation_info.R script")
}

# Save output
message("Saving outputs")
write_csv(df, file = out_csv)

# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)
