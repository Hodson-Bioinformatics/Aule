#!/usr/bin/env Rscript
#
# Usage:
#   get_bams.R <input_table> <output_table> [<library_id_column>]
#
# Notes:
#   - <library_id_column> can be a name or the number 1. If it is the
#     number 1, the script assumes that the input file has no header
#     and only has a single column. The default value is 'library_id'
#     (i.e., the script assumes a file header).
#
#   - The column names are prefixed with 'al.', 'lc.', and 'lb.' for
#     designating alignment-, libcore-, and library-specific columns.
#     It also helps with ensuring no collisions when joining the
#     obtained information with the input table.
#


# Load packages -----------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages({
  
  message("Loading packages...")
  
  required_packages <- c(
   "tidyverse"
  )
  
  installed_packages <- rownames(installed.packages())
  
  for (pkg in required_packages){
    if (!pkg %in% installed_packages) {
      message(paste0("Installing ", pkg, "..."))
      install.packages(pkg, repos = "https://mirror.rcg.sfu.ca/mirror/CRAN/")
    }
    library(pkg, character.only = TRUE)
  }
}))

# Determine arguments -----------------------------------------------------

# Parse command-line arguments (when running non-interactively)
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("input_table", "output_table", "library_id_column")
args <- setNames(args, arg_names[1:length(args)])

# Set default value
if (is.null(args$library_id_column)) args$library_id_column <- "library_id"


# Check if table or single-column list
is_single_column <- args$library_id_column == "1"



# Load data ---------------------------------------------------------------

message("Loading data...")

if (is_single_column) {
  libraries <- read_tsv(args$input_table, col_names = "library_id",
                        col_types = cols())
  library_id_column = "library_id"
} else {
  libraries <- read_tsv(args$input_table, col_types = cols())
  library_id_column = args$library_id_column
}

# Print error if args$library_id_column isn't one of the colnames
if(!library_id_column %in% colnames(libraries)) stop(paste0('There is no column called "', library_id_column, '" in the input table. \n Check usage.'))
  

# Load get_bams() function--------------------------------------------------

this_file <- gsub("--file=", "", commandArgs()[grepl("--file", commandArgs())])
get_function <- gsub(".R", "_function.R", this_file)
source(get_function)

output_table <- get_bams(libraries, library_id_column)

message("Writing final table to file")

write_tsv(output_table, args$output_table)


