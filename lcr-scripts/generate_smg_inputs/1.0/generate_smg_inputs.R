### Generate inputs for SMG tools ###
# Uses file defining sample set to generate maf
# and prepare it to be used with SMG tools.

#!/usr/bin/env Rscript
#
# Usage:
#   Rscript generate_smg_inputs.R <path/to/master/maf> <path/to/sample_sets> <path/to/output/folder> <case_set> <mode> <include non-coding>
#
# Notes:
#   This script is intended for use with the SMG modules in LCR-modules (MutSig2CV, dNdS, HotMAPS,
#   OncodriveFML).
#   It expects to be provided with the tab-deliminated file where sample subsets
#   for a particular analysis are specified, where first column (sample_id/Tumor_Sample_Barcode)
#   defines the unique sample ID, and each column indicates whether this ID is included (1) or not (0)
#   in a particular subset. The column name for the subset will be used as the naming of the
#   output maf file at the user-provided location.
#

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressWarnings(
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})
)

# Determine arguments -----------------------------------------------------

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("master_maf", "all_sample_sets", "output_path", "case_set", "mode", "include_non_coding")
# if there are multiple maf files passed, collapse them into one list
args = c(
        list(unlist(args[1:(length(args)-5)])),
        args[(length(args)-4):length(args)]
        )
args <- setNames(args, arg_names[1:length(args)])
args$master_maf = as.list(args$master_maf)

# Print args for debugging
print(paste("master_maf:",args$master_maf,
            "all_sample_sets:",args$all_sample_sets,
            "output_path:",args$output_path,
            "case_set",args$case_set,
            "mode",args$mode,
            "include_non_coding",args$include_non_coding))

# Ensure consistent naming of the sample ID column-------------------
if (file.exists(args$all_sample_sets)) {
  full_case_set = suppressMessages(read_tsv(args$all_sample_sets))
} else {
  message(paste("Warning: case set is requested, but the case set file", full_case_set_path, "is not found."))
  stop("Exiting because did not found case set-defining file")
}

full_case_set =
  full_case_set %>% rename_at(vars(matches(
    "sample_id", ignore.case = TRUE
  )),
  ~ "Tumor_Sample_Barcode")

# get case set as defined in the file
this_subset_samples =
  full_case_set %>%
  dplyr::filter(!!sym(args$case_set) == 1) %>%
  pull(Tumor_Sample_Barcode)


# Load master maf and get mutations for the maf subset-------------------
message("Loading master maf and finding available data for samples in requested subset...")
if (length(args$master_maf)>1){
  message("More than one maf file is supplied. Concatenating them into single file.")
  master_maf =
    tibble(filename = args$master_maf) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename, # read files into
                             ~ read_tsv(args$master_maf, col_types = cols())) # a new data column
    ) %>%
    unnest(cols = c(file_contents)) %>%
    select(-filename)
} else {
  master_maf = suppressWarnings(read_tsv(args$master_maf, col_types = cols()))
}

subset_maf =
  master_maf %>%
  dplyr::filter(Tumor_Sample_Barcode %in% this_subset_samples)

# subset for coding only if user requested
if (args$include_non_coding) {
  message("Proceeding with non-coding mutations...")
} else{
  message("Excluding non-coding mutations...")
  coding_class = c("Frame_Shift_Del",
                "Frame_Shift_Ins",
                "In_Frame_Del",
                "In_Frame_Ins",
                "Missense_Mutation",
                "Nonsense_Mutation",
                "Nonstop_Mutation",
                "Silent",
                "Splice_Region",
                "Splice_Site",
                "Targeted_Region",
                "Translation_Start_Site")
    subset_maf =
      subset_maf %>%
      dplyr::filter(Variant_Classification %in% coding_class)
}

# report any missing samples
missing_samples = setdiff(this_subset_samples,
                          unique(subset_maf$Tumor_Sample_Barcode))

if (length(missing_samples)==0) {
  message(paste("Success! Found mutations for all samples.", length(this_subset_samples), "patients will be used in the analysis"))
} else {
  message(paste("WARNING:", length(missing_samples), "will not be available for the analysis."))
  message("Did not find mutations for these samples in the master input maf:")
  message(missing_samples)
}

# Format maf according to the requirements of each individual tool --------------------------------------
message(paste("preparing maf file to be used with", args$mode))
# MutSig2CV
if (args$mode == "MutSig2CV") {
  # check that the maf file is in grch37-based coordinates
  if (grepl("38", subset_maf$NCBI_Build[1])) {
    message("Requested mode is MutSig2CV, but the supplied file is in the hg38-based coordinates.")
    message("Unfortunatelly, MutSig only works for grch37-based maf files.")
    stop("Please supply the mutation data in grch37-based version.")
  }

  subset_maf =
    subset_maf %>%
    rename("chr"="Chromosome",
           "pos"="Start_Position",
           "gene"="Hugo_Symbol",
           "patient"="Tumor_Sample_Barcode",
           "ref_allele"="Reference_Allele",
           "newbase"="Tumor_Seq_Allele2",
           "type"="Variant_Classification",
           "classification"="Variant_Type")

  subset_maf =
  subset_maf %>%
    select(chr, pos, gene, patient, ref_allele, newbase, type, classification)

  grouping_column = "patient"

}

# dNdS
if (args$mode == "dNdS") {
  # check that the maf file is in grch37-based coordinates
  if (grepl("38", subset_maf$NCBI_Build[1])) {
    message("Requested mode is dNdS, but the supplied file is in the hg38-based coordinates.")
    message("Unfortunatelly, dNdS is configured to only work for grch37-based maf files.")
    stop("Please supply the mutation data in grch37-based version.")
  }

  subset_maf =
    subset_maf %>%
      select(Tumor_Sample_Barcode,
              Chromosome,
              Start_Position,
              Reference_Allele,
              Tumor_Seq_Allele2) %>%
      `names<-`(c("sampleID",
                  "chr",
                  "pos",
                  "ref",
                  "mut"))

  grouping_column = "sampleID"

}

# prepare maf file contents for documentation purposes
contents = subset_maf %>%
  group_by(across(all_of(grouping_column))) %>%
  summarise(N_mutations=n()) %>%
  ungroup %>%
  mutate(non_coding_included=args$include_non_coding)

# Output --------------------------------------------------------------------------------------

message("Writing resulting maf to file...")
write_tsv(subset_maf, paste0(args$output_path, "/", args$case_set, ".maf"))
write_tsv(contents, paste0(args$output_path, "/", args$case_set, ".maf.content"))

message("DONE!")
