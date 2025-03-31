# Author: Natasha Cutmore
# Contributor:  Ilias Moutsopoulos
# Date: 09/12/2022
# Revised: 30/08/2024
# Title: SVs output
# Description: Merge GRIDSS and Manta output and make output data table

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
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GenomicRanges))

# 1.B Parse Snakemake object
merged_csv  <- snakemake@input[["merged_csv"]]
ROI_with_IG <- snakemake@input[["roi_annotations"]]
outfile     <- snakemake@output[["csv"]]
log_stdout <- snakemake@log[["stdout"]]
log_stderr <- snakemake@log[["stderr"]]

# Capture output
file_stderr <- file(log_stderr, open = "at")
sink(file_stderr, type = "message", append = TRUE)
sink(log_stdout, type = "output", append = TRUE)

# 1.C Load Functions:
`%!in%` <- Negate(`%in%`)

blank_tab <- tibble(
                    ID = NA_character_,
                    caller.ID = NA_character_,
                    Partner1_chr = NA_character_, 
                    Partner1_pos = NA_integer_, 
                    Partner1_gene = NA_character_,
                    Partner2_chr = NA_character_, 
                    Partner2_pos = NA_integer_, 
                    Partner2_gene = NA_character_,
                    VAF = NA_integer_,
                    RP = NA_integer_,
                    SR = NA_integer_,
                    Found_in_GRIDSS = NA_integer_,
                    Found_in_Manta = NA_integer_ 
                )

tibble_2_GRanges <- function(chrom, start, end, names = NA){
    
    GRanges(seqnames= c(chrom),
            IRanges(start = as.numeric(c(start)),
                    end = as.numeric(c(end))),
            names = as.character(c(names)))
    
}

format_SVs <- function(csv.tab){
    
    if(!("VAF" %in% colnames(csv.tab))) csv.tab["VAF"] <- NA_integer_
    
    Partner1.gr <- tibble_2_GRanges(csv.tab$Partner1_chr, 
                                    csv.tab$Partner1_pos, 
                                    csv.tab$Partner1_pos)
    
    Partner2.gr <- tibble_2_GRanges(csv.tab$Partner2_chr, 
                                    csv.tab$Partner2_pos, 
                                    csv.tab$Partner2_pos)
    
    Partner1.ROI.olaps <- findOverlaps(Partner1.gr, ROI.annotations.gr, maxgap = 100)
    
    Partner2.ROI.olaps <- findOverlaps(Partner2.gr, ROI.annotations.gr, maxgap = 100)
    
    p1_qh <- queryHits(Partner1.ROI.olaps)
    p1_sh <- subjectHits(Partner1.ROI.olaps)
    p2_qh <- queryHits(Partner2.ROI.olaps)
    p2_sh <- subjectHits(Partner2.ROI.olaps)
    
    annotated.tab <- csv.tab %>% 
        tibble::rowid_to_column(var = "rowid") %>%
        dplyr::mutate(Partner1_gene = ifelse(rowid %in% p1_qh, ROI.annotations.gr$names[p1_sh[match(rowid, p1_qh)]], "Partner unknown"),
                      Partner2_gene = ifelse(rowid %in% p2_qh, ROI.annotations.gr$names[p2_sh[match(rowid, p2_qh)]], "Partner unknown")) %>%
        dplyr::distinct() %>%
        dplyr::select(ID, caller.ID, Partner1_chr, Partner1_pos, Partner1_gene, Partner2_chr, Partner2_pos, Partner2_gene, VAF, SR, RP) %>%
        dplyr::mutate(Partner1_pos = as.numeric(Partner1_pos),
                      Partner2_pos = as.numeric(Partner2_pos),
                      PartnerA_chr = dplyr::case_when(Partner1_chr > Partner2_chr ~ Partner1_chr,
                                                      Partner1_chr == Partner2_chr ~ Partner1_chr,
                                                      TRUE ~ Partner2_chr),
                      PartnerA_pos = dplyr::case_when(Partner1_chr > Partner2_chr ~ Partner1_pos,
                                                      Partner1_chr == Partner2_chr ~ Partner1_pos,
                                                      TRUE ~ Partner2_pos),
                      PartnerA_gene = dplyr::case_when(Partner1_chr > Partner2_chr ~ Partner1_gene,
                                                       Partner1_chr == Partner2_chr ~ Partner1_gene,
                                                       TRUE ~ Partner2_gene),
                      PartnerB_chr = dplyr::case_when(Partner1_chr > Partner2_chr ~ Partner2_chr,
                                                      Partner1_chr == Partner2_chr ~ Partner2_chr,
                                                      TRUE ~ Partner1_chr),
                      PartnerB_pos = dplyr::case_when(Partner1_chr > Partner2_chr ~ Partner2_pos,
                                                      Partner1_chr == Partner2_chr ~ Partner2_pos,
                                                      TRUE ~ Partner1_pos),
                      PartnerB_gene = dplyr::case_when(Partner1_chr > Partner2_chr ~ Partner2_gene,
                                                       Partner1_chr == Partner2_chr ~ Partner2_gene,
                                                       TRUE ~ Partner1_gene)) %>%
        dplyr::select(ID, caller.ID, PartnerA_chr, PartnerA_pos, PartnerA_gene, PartnerB_chr, PartnerB_pos, PartnerB_gene, VAF, SR, RP) %>%
        dplyr::rename(Partner1_chr = PartnerA_chr,
                      Partner1_pos = PartnerA_pos,
                      Partner1_gene = PartnerA_gene,
                      Partner2_chr = PartnerB_chr,
                      Partner2_pos = PartnerB_pos,
                      Partner2_gene = PartnerB_gene) %>%
        dplyr::filter(Partner1_chr %in% ROI.annotations$chromosome_name |
                      Partner2_chr %in% ROI.annotations$chromosome_name)
    
    
}

annotate_SVs <- function(input.csv){
    
    csv.in <- readr::read_csv(input.csv, show_col_types = FALSE) %>%
        dplyr::filter(!is.na(Partner1_chr))
    
    if(nrow(csv.in) == 0){
        
        filename <- basename(input.csv)
        Normal <- unlist(strsplit(filename, "--"))[2]
        Tumour <- unlist(strsplit(filename, "--"))[1]
        
        annotated.tab <- blank_tab %>%
            dplyr::mutate(Indiv = Tumour) %>%
            dplyr::select(Indiv, everything())
        
    } else {
        
        filename <- basename(input.csv)
        Normal <- unlist(strsplit(filename, "--"))[2]
        Tumour <- unlist(strsplit(filename, "--"))[1]
        
        csv.in <- read_csv(input.csv, show_col_types = FALSE)
        
        annotated.tab <- tibble(format_SVs(csv.in)) %>%
            dplyr::mutate(Indiv = Tumour) %>%
            dplyr::select(Indiv, everything())
        
    }
    
    return(annotated.tab)
    
}


# 3. Wrangle

# 3.A Get Region of Interest info
ROI.annotations <- read.csv(ROI_with_IG) %>% 
    dplyr::rename(chromosome_name = chr)
ROI.annotations$names <- as.character(ROI.annotations$names)
ROI.annotations$chromosome_name <- as.character(ROI.annotations$chromosome_name)

ROI.annotations.gr <- tibble_2_GRanges(ROI.annotations$chromosome_name, 
                                       ROI.annotations$start, 
                                       ROI.annotations$end,
                                       ROI.annotations$names)

# 3.B Annotate, filter and tidy
annotatedSVs <- annotate_SVs(merged_csv)

write_csv(annotatedSVs, outfile) 

# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)
