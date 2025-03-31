# title: "01_lymphgen_SAF"
# author: "Natasha Cutmore"
# date: "31/03/2023"
# output: tab delim .txt
# Aim: Generate LymphGen sample annotation file.

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
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))

# Load Options
option_list = list(
  make_option(c("--by_sv_samples"), type="character", default=NULL,
              help="by_sv_sample.csv files, .csv", metavar="character"),
  make_option(c("--samples"), type="character", default=NULL,
              help="by_sample.csv files, .csv", metavar="character"),
  make_option(c("--out_txt"), type="character", default=NULL,
              help="Output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

samples       <- opt$samples
samples_svs   <- opt$by_sv_samples
out_txt       <- opt$out_txt

# Load Data
samples.csv <- read.csv(samples) %>% dplyr::select(sample_id)
samples_svs.csv <- read.csv(samples_svs)

# Wrangle into LG sample annotation file
output.table <- samples.csv %>%
  dplyr::left_join(samples_svs.csv, by = "sample_id") %>%
  dplyr::rename(Sample.ID = sample_id) %>%
  rowwise() %>%
  dplyr::mutate(Copy.Number = 0,
                BCL2.transloc = case_when(BCL2 == 1 ~ 1, TRUE ~0), 
                BCL6.transloc = case_when(BCL6 == 1 ~1, TRUE ~0)  
                ) %>%
  ungroup() %>%
  distinct() %>%
  na.omit() %>%
  dplyr::select(Sample.ID, Copy.Number, BCL2.transloc, BCL6.transloc)


write_delim(output.table, out_txt, delim = "\t") 

