# Author: Ilias Moutsopoulos
# Date: 02/03/2023
# Title: Join data
# Description: Create table of information included in a run

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
suppressPackageStartupMessages(library(optparse))

# Load Arguments
option_list = list(
  make_option(c("--samples"), type="character", default=NULL, 
              help="Sample table for run", metavar="character"),
  make_option(c("--in_csv"), type="character", default=NULL, 
              help="Path of csv with extra information", metavar="character"),
  make_option(c("--out_csv"), type="character", default=NULL, 
              help="Output csv path", metavar="character"),
  make_option(c("--column_name"), type="character", default=NULL, 
              help="Name of column to extract data by, must be present in both inputs", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

samples     <- opt$samples
in_csv      <- opt$in_csv
out_csv     <- opt$out_csv
column_name <- opt$column_name

# Create, join, save table
df <- read.table(samples, sep = "\t", header = TRUE)

if (length(unique(df[[column_name]])) != nrow(df)) 
  df <- df %>% distinct(!!sym(column_name))

df %>%
  left_join(read.csv(in_csv), by = column_name) %>%
  write.csv(out_csv, row.names = FALSE)
