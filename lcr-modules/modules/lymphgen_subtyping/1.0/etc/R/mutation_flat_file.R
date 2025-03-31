# title: "01_HMRN_format_matrix"
# author: "Natasha Cutmore"
# date: "31/03/2023"
# output: html_document
# Aim: Generate HMRN classifier-ready matrix, for subtyping.

## If R is installed with conda, this will configure the 
## library paths to only point to the conda R package directory
local({
  conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
  rhome <- normalizePath(R.home())
  if (!is.na(conda) && startsWith(rhome, conda)) {
    .libPaths(R.home("library"), include.site = FALSE)
  }
})

# Load Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))

# Load Options
option_list = list(
  make_option(c("--samples"), type="character", default=NULL,
              help="by_sample.csv files, .csv", metavar="character"),
  make_option(c("--mutation_file"), type="character", default=NULL,
              help="by_mutation.csv files, .csv", metavar="character"),
  make_option(c("--out_mff"), type="character", default=NULL,
              help="Output mutation flat file, .txt", metavar="character"),
  make_option(c("--out_mgl"), type="character", default=NULL,
              help="Output mutation gene list, .txt", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mutation_file <-  opt$mutation_file
out_mff       <-  opt$out_mff
out_mgl       <-  opt$out_mgl

# Load Data
input.csv <- read.csv(mutation_file) 

# Get input mutations and wrangle into LG mutation flat file *from ckrushton, LGenIC https://github.com/ckrushton/LGenIC
`%!in%` <- Negate(`%in%`)

silent.variants <- c("Silent", "5'UTR", "5'Flank", "Intron", "3'UTR")
trunc.variants <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")
missense.variants <- c("In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
special.genes <-c("NOTCH1", "NOTCH2", "CD79B", "EZH2")

input.processed.df <- input.csv %>% 
  dplyr::filter(!grepl("^IG", Hugo_Symbol)) %>%
  dplyr::rename(Sample = sample_id,
                ENTREZ.ID = Entrez_Gene_Id) %>%
  dplyr::mutate(Type = case_when(
                    (Hugo_Symbol == "MYD88" & grepl("p.Leu273Pro", HGVSp)) ~ "L265P",
                    (Variant_Classification %in% silent.variants) ~ "Synon",
                    (Variant_Classification %in% trunc.variants) ~ "TRUNC",
                    (Variant_Classification %in% missense.variants) ~ "MUTATION")) %>%
  dplyr::mutate(ENTREZ.ID = ifelse(
                    (Hugo_Symbol == "RP11-211G3.2" & is.na(ENTREZ.ID)), 604, ENTREZ.ID),
                Keep_mutation = case_when(
                    (Hugo_Symbol == "NOTCH1" & CHROM == "chr9" & between(POS, 136494433, 136497003) & Type == "TRUNC") ~1,
                    (Hugo_Symbol == "NOTCH2" & CHROM == "chr1" & between(POS, 119911553, 119916527) & Type == "TRUNC") ~1,
                    (Hugo_Symbol == "CD79B" & CHROM == "chr17" & between(POS, 63929439, 63929812) & Type == "TRUNC") ~1,
                    (Hugo_Symbol == "CD79B" & CHROM == "chr17" & between(POS, 63928740, 63929440) & Type == "MUTATION") ~1,
                    (Hugo_Symbol == "EZH2" & CHROM == "chr7" & between(POS, 148809146, 148811672)) ~1,
                    (Hugo_Symbol %!in% special.genes) ~1,
                    TRUE ~ 0)) %>%
  dplyr::filter(!grepl("^intron", Consequence),
                !grepl("^intergenic", Consequence),
                !is.na(Type),
                Keep_mutation == 1) %>%
  dplyr::select(Sample, ENTREZ.ID, 
                 Type) %>%
  distinct()


write_delim(input.processed.df, out_mff, delim = "\t")


# Wrangle into LG mutation gene list 
mut_gene_list <- input.processed.df %>%
  dplyr::select(ENTREZ.ID) %>%
  dplyr::distinct() %>%
  na.omit()

write_delim(mut_gene_list, out_mgl, delim = "\t")

#####Sandbox -for trouble shooting
# save.image("/servers/bigscratch/hodson/DIRECT/LCR/main/runs/20230314_NHC_testrun/.RDdata")
# load("/servers/bigscratch/hodson/DIRECT/LCR/main/runs/20230314_NHC_testrun/.RDdata") ## Loads as environment
# open terminal, conda activate env from error, then execute above.
