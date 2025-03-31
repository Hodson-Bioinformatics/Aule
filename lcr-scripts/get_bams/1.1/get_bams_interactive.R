#!/usr/bin/env Rscript
#
# Usage:
#   get_bams_interactive.R "<library_id>"
#
# Notes:
#   - <library_id> may be a single GSC library ID or a comma-separated list
#     of IDs, e.g. get_bams_interactive.R "D85283,D98729". 
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
# Error if no library_id is specified. 
if (length(args) != 1) stop("Please specify a single library ID or comma-separated list of library IDs. \n\nUsage: get_bams_interactive.R \"D85283,D98729\"")
arg_names <- c("library_id")
args <- setNames(args, arg_names[1:length(args)])



# Load data ---------------------------------------------------------------

message("Loading data...")

libraries <- data.frame(library_id = unlist(str_split(str_remove_all(args$library_id, " "), pattern = ",")))


# Load get_bams() function--------------------------------------------------

this_file <- gsub("--file=", "", commandArgs()[grepl("--file", commandArgs())])
get_function <- gsub("_interactive.R", "_function.R", this_file)
source(get_function)

output_table <- get_bams(libraries, "library_id") 

cat(format_tsv(select(output_table, library_id, al.merge_status, al.data_path, lb.gsc_external_id, lb.protocol_name)))



