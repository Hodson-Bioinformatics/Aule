# Author: Natasha Cutmore
# Contributor:  Ilias Moutsopoulos
# Date: 02/09/2024
# Title: merge caller output
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
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(StructuralVariantAnnotation))

# 1.B Parse Snakemake object
manta_vcf  <- snakemake@input[["manta_vcf"]]
gridss_vcf <- snakemake@input[["gridss_vcf"]]
out_csv    <- snakemake@output[["csv"]]
log_stdout <- snakemake@log[["stdout"]]
log_stderr <- snakemake@log[["stderr"]]

# Capture output
file_stderr <- file(log_stderr, open = "at")
sink(file_stderr, type = "message", append = TRUE)
sink(log_stdout, type = "output", append = TRUE)

# 1.C Load Functions
`%!in%` <- Negate(`%in%`)

blank_tab <- tibble(
    ID = NA_character_,
    caller.ID = NA_character_,
    Partner1_chr = NA_character_, 
    Partner1_pos = NA_integer_,
    REF = NA_character_,
    ALT = NA_character_,
    Partner2_chr = NA_character_, 
    Partner2_pos = NA_integer_,
    QUAL = NA_integer_,
    FILTER = NA_character_,
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

gridss_vcf_2_tibble <- function(input.vcf){
    
    vcf.in <- vcfR::read.vcfR(input.vcf) 
    
    if(nrow(vcf.in) == 0){
        
        filename <- basename(input.vcf)
        Normal <- unlist(strsplit(filename, "--"))[2]
        Tumour <- unlist(strsplit(filename, "--"))[1]
        
        vcf.tab <- blank_tab %>% 
            dplyr::mutate(Indiv = Tumour) %>%
            dplyr::select(Indiv, everything())
        
        
    } else {
        
        vcf.in <- vcfR::read.vcfR(input.vcf) %>%
            vcfR::vcfR2tidy()
        
        filename <- basename(input.vcf)
        Normal <- unlist(strsplit(filename, "--"))[2]
        Tumour <- unlist(strsplit(filename, "--"))[1]
        
        vcf.tab <- tibble(vcf.in$gt) %>%
            dplyr::filter(grepl(Tumour, Indiv)) %>%
            dplyr::select((-c("POS", "ChromKey"))) %>%
            cbind(vcf.in$fix) %>%
            dplyr::rename(VAF = TAF) %>%
            dplyr::select(-c(starts_with("gt_"))) %>%
            tidyr::separate(ALT, 
                            sep = "[\\[\\]:]+",
                            into = c("temp1", "Partner2_chr", "Partner2_pos", "temp2"),
                            remove = FALSE) %>%
            dplyr::filter(!is.na(temp1), 
                          !is.na(Partner2_chr), 
                          !is.na(Partner2_pos), 
                          !is.na(temp2)
            ) %>%
            dplyr::mutate(caller.ID = NA_character_,
                          RP = as.numeric(RP),
                          SR = as.numeric(SR)) %>%
            dplyr::select(Indiv, ID, caller.ID, CHROM, POS, REF...4, ALT, Partner2_chr, Partner2_pos, 
                          QUAL, FILTER, SR, RP, VAF) %>%
            dplyr::rename(Partner1_chr = CHROM,
                          Partner1_pos = POS,
                          REF = REF...4) %>%
            dplyr::mutate(Found_in_GRIDSS = 1,
                          Found_in_Manta = 0)
        
    }
    return(vcf.tab)
    
    
}

manta_vcf_2_tibble <- function(input.vcf){
    
    vcf.in <- vcfR::read.vcfR(input.vcf) 
    
    if(nrow(vcf.in) == 0){
        
        
        filename <- basename(input.vcf)
        Normal <- unlist(strsplit(filename, "--"))[2]
        Tumour <- unlist(strsplit(filename, "--"))[1]
        
        vcf.tab <- blank_tab %>% 
            dplyr::mutate(Indiv = Tumour) %>%
            dplyr::select(Indiv, everything())
        
        
    } else {
        
        vcf.in <- vcfR::read.vcfR(input.vcf) %>%
            vcfR::vcfR2tidy()
        
        filename <- basename(input.vcf)
        Normal <- unlist(strsplit(filename, "--"))[2]
        Tumour <- unlist(strsplit(filename, "--"))[1]
        
        
        vcf.tab <- tibble(vcf.in$gt) %>%
            dplyr::filter(grepl(Tumour, Indiv)) %>%
            dplyr::select((-c("POS", "ChromKey"))) %>%
            cbind(vcf.in$fix) %>%
            dplyr::rename(VAF = gt_VAF,
                          RP = gt_PR,
                          SR = gt_SR) %>%
            dplyr::select(-c(starts_with("gt_"))) %>%
            tidyr::separate(ALT, into = c("ALT_base",
                                          "Partner2_chr",
                                          "Partner2_pos"),
                            remove = FALSE) %>%
            tidyr::separate(SR, sep = "[,]+", into = c("SR_REF", "SR"), remove = TRUE) %>%
            tidyr::separate(RP, sep = "[,]+", into = c("RP_REF", "RP"), remove = TRUE) %>%
            dplyr::filter(SVTYPE =="BND") %>%
            dplyr::mutate(caller.ID = NA_character_,
                          RP = as.numeric(RP),
                          SR = as.numeric(SR)) %>%
            dplyr::select(Indiv, ID, caller.ID, CHROM, POS, REF, ALT, Partner2_chr, Partner2_pos, 
                          QUAL, FILTER, SR, RP, VAF) %>%
            dplyr::rename(Partner1_chr = CHROM,
                          Partner1_pos = POS) %>%
            dplyr::mutate(Found_in_GRIDSS = 0,
                          Found_in_Manta = 1)
        
    }
    return(vcf.tab) 
    
}

olap_svs <- function(GRIDSS_infile, Manta_infile){
    
    # read in vcf files 
    gridss.vcf <- VariantAnnotation::readVcf(GRIDSS_infile)
    manta.vcf  <- VariantAnnotation::readVcf(Manta_infile)
    
    # convert to breakpoint ranges objects
    gridss.gr <- breakpointRanges(gridss.vcf)
    manta.gr  <- breakpointRanges(manta.vcf)
    
    # find overlaps
    olaps <- findBreakpointOverlaps(gridss.gr, manta.gr,
                                    maxgap = 1000,
                                    sizemargin = 0.25,
                                    restrictMarginToSizeMultiple = 0.5)
    
    qh <- queryHits(olaps)
    sh <- subjectHits(olaps)
    
    # retrieve overlap data
    gridss.tab <- gridss_vcf_2_tibble(GRIDSS_infile) %>% 
        tibble::rowid_to_column(var = "rowid") %>%
        dplyr::mutate(Found_in_both = ifelse(rowid %in% qh, 
                                             manta.gr@ranges@NAMES[sh[match(rowid, qh)]], 
                                             NA_character_)) %>%
        dplyr::filter(grepl("o$", ID)) # GRIDSS has every SV represented twice
    
    # tidy and format into output table
    manta.tab  <- manta_vcf_2_tibble(Manta_infile) %>%
        dplyr::filter(ID %!in% gridss.tab$Found_in_both)
    
    all <- gridss.tab %>%
        dplyr::mutate(Found_in_Manta = ifelse(!is.na(Found_in_both), 1, 0)) %>%
        dplyr::select(-c(Found_in_both, rowid)) %>%
        rbind(manta.tab)
    
    return(all)
    
}

# 3. Wrangle
output <- olap_svs(GRIDSS_infile = gridss_vcf, Manta_infile = manta_vcf)

write.csv(output, out_csv, row.names = FALSE, quote = FALSE)


# Close log file connections
sink(type = "output")
sink(type = "message")
close(file_stderr)
