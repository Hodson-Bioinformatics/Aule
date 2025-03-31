# Author: Natasha Cutmore
# Contributor:
# Date: 30/08/2024
# Title: Summarise SVs output
# Description: Summarise output for sv_merge into per-sample table format

## If R is installed with conda, this will configure the 
## library paths to only point to the conda R package directory
local({
    conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
    rhome <- normalizePath(R.home())
    if (!is.na(conda) && startsWith(rhome, conda)) {
        .libPaths(R.home("library"), include.site = FALSE)
    }
})

# 1.A Load Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

# 1.B Parse snakemake object
csvs       <- snakemake@input[["csvs"]]
samples    <- snakemake@input[["samples"]]
SR_filter  <- snakemake@params[["SR_filter"]]
RP_filter  <- snakemake@params[["RP_filter"]]
out_csv    <- snakemake@output[["csv"]]
log_stdout <- snakemake@log[["stdout"]]
log_stderr <- snakemake@log[["stderr"]]

# Capture output
file_stderr <- file(log_stderr, open = "at")
sink(file_stderr, type = "message", append = TRUE)
sink(log_stdout, type = "output", append = TRUE)


# 1.C Load Functions
`%!in%` <- Negate(`%in%`)

# 2. Load data
samples <- readr::read_delim(samples, show_col_types = FALSE) %>%
    dplyr::filter(direct_type != "Germline")

# 3. Summarise SVs
# Read, tidy and combine all mafs
csv_df <- tibble()
for(csv in csvs){
    temp <- readr::read_csv(csv, show_col_types = FALSE) %>%
        dplyr::mutate(sample_id = unlist(strsplit(basename(csv), "--"))[1], .before = 1) %>%
        dplyr::select(sample_id, starts_with("Partner"), SR, RP) %>%
        dplyr::filter(!(grepl("^Ig", Partner1_gene) & grepl("Ig", Partner2_gene)),
                      !(grepl("^Partner", Partner1_gene) & grepl("^Ig", Partner2_gene)),
                      !(grepl("^Partner", Partner1_gene) & grepl("^Partner", Partner2_gene)),
                      !(grepl("^Ig", Partner1_gene) & grepl("^Partner", Partner2_gene)),
                      !(grepl("BCL2", Partner1_gene) & grepl("BCL2", Partner2_gene)),
                      (SR >= SR_filter & RP >0) | (SR > 0 & RP >= RP_filter)) %>%
        dplyr::group_by(sample_id, Partner1_gene, Partner2_gene) %>%
        dplyr::summarise(N = n()) %>%
        dplyr::mutate(N = ifelse(N >= 1, 1, 0)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = c("Partner1_gene", "Partner2_gene"),
                           values_from = "N",
                           names_glue = "{Partner1_gene}_{Partner2_gene}") %>%
        unnest() %>%
        distinct() %>%
        dplyr::mutate(across(where(is.numeric), ~replace_na(.,0)))
        
    if (nrow(csv_df) == 0) {
        csv_df <- temp
    } else {
        csv_df <- csv_df %>%
            full_join(temp)
    }
}

all_samples <- expand.grid(sample_id = unique(samples$sample_id))

svs_summarised <- tibble(sample_id = NA_character_,
                         BCL2_IgH = NA_integer_,
                         BCL2_IgK = NA_integer_,
                         BCL2_IgL = NA_integer_,
                         `BCL2_Partner unknown` = NA_integer_,
                         BCL6_IgH = NA_integer_,
                         BCL6_IgK = NA_integer_,
                         BCL6_IgL = NA_integer_,
                         `BCL6_Partner unknown` = NA_integer_,
                         BCL6_BCL6 = NA_integer_,
                         MYC_IgH = NA_integer_,
                         MYC_IgK = NA_integer_,
                         MYC_IgL = NA_integer_,
                         `MYC_Partner unknown` = NA_integer_,
                         MYC_MYC = NA_integer_) %>%
    dplyr::right_join(csv_df) 

full_df <- all_samples %>%
    dplyr::full_join(svs_summarised)
    

write_csv(full_df, out_csv)

# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)
