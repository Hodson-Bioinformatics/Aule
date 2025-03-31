# Author: Natasha Cutmore
# Contributor:  Ilias Moutsopoulos
# Date: 09/12/2022
# Title: SVs output
# Description: Take IgCaller and intersected GRIDSS/Manta output and make output data table

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
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(optparse))

# 1.B Load Arguments
option_list = list(
  make_option(c("-i", "--IgCaller_infile"), type="character", default=NULL,
              help="IgCaller *_filtered.tsv", metavar="character"),
  make_option(c("-s", "--SURVIVOR_infile"), type="character", default=NULL,
              help="Merged SURVIVOR output file, .vcf", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="Output merged csv", metavar="character"),
  make_option(c("-g", "--GRIDSS_infile"), type="character", default=NULL,
              help="GRIDSS filtered .vcf", metavar="character"),
  make_option(c("-m", "--Manta_infile"), type="character", default=NULL,
              help="Manta filtered .vcf", metavar="character"),
  make_option(c("-r", "--ROI_with_IG"), type="character", default=NULL,
              help="Regions of interest csv including IG loci", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# 1.C Load Functions:

tibble_2_GRanges <- function(chrom, start, end, names = NA){
  
  GRanges(seqnames= c(chrom),
          IRanges(start = as.numeric(c(start)),
                  end = as.numeric(c(end))),
          names = as.character(c(names)))
  
}

adjust_overlaps <- function(partner.gr, olaps.gr){
  qh <- queryHits(olaps.gr)
  sh <- subjectHits(olaps.gr)
  missing <- seq_len(length(partner.gr) - length(olaps.gr))
  missing_query_1 <- setdiff(seq_along(partner.gr), qh)
  qh <- c(qh, missing)
  sh <- c(sh, rep(NA, length(missing)))
  list(qh = qh, sh = sh)
}

format_SVs <- function(vcf.tab){
  
  if(!("VAF" %in% colnames(vcf.tab))) vcf.tab["VAF"] <- NA_integer_
  
  Partner1.gr <- tibble_2_GRanges(vcf.tab$CHROM, 
                                  vcf.tab$POS, 
                                  vcf.tab$POS)
  
  Partner2.gr <- tibble_2_GRanges(vcf.tab$Partner2_chr, 
                                  vcf.tab$Partner2_pos, 
                                  vcf.tab$Partner2_pos)
  
  Partner1.ROI.olaps <- findOverlaps(Partner1.gr, ROI.annotations.gr)
  
  Partner2.ROI.olaps <- findOverlaps(Partner2.gr, ROI.annotations.gr)
  
  qh_sh_1 <- adjust_overlaps(Partner1.gr, Partner1.ROI.olaps)
  P1.ROI.qh <- qh_sh_1$qh
  P1.ROI.sh <- qh_sh_1$sh
  
  qh_sh_2 <- adjust_overlaps(Partner2.gr, Partner2.ROI.olaps)
  P2.ROI.qh <- qh_sh_2$qh
  P2.ROI.sh <- qh_sh_2$sh
  
  annotated.tab <- vcf.tab[P1.ROI.qh,] %>%
    dplyr::mutate(Partner1_gene = ROI.annotations.gr$names[P1.ROI.sh],
                  Partner2_gene = ROI.annotations.gr$names[P2.ROI.sh]) %>%
    dplyr::distinct() %>%
    dplyr::select(ID, CHROM, POS, Partner1_gene, Partner2_chr, Partner2_pos, Partner2_gene, VAF) %>%
    dplyr::rename(Partner1_chr = CHROM,
                  Partner1_pos = POS) %>%
    dplyr::mutate(SV_ID = dplyr::case_when(Partner1_chr > Partner2_chr ~ substr(paste(Partner1_chr, Partner1_pos, sep = "_"), 1, 13),
                                           Partner2_chr > Partner1_chr ~ substr(paste(Partner2_chr, Partner2_pos, sep = "_"), 1, 13),
                                           Partner1_chr == Partner2_chr ~substr(paste(Partner1_chr, Partner1_pos, sep = "_"), 1, 13))) 
  
  
}

survivor_vcf_2_tibble <- function(input.vcf){
  
  vcf.in <- vcfR::read.vcfR(input.vcf) 
  
  if(nrow(vcf.in) == 0){
    
    annotated.tab <- tibble(
      SURVIVOR.Partner1_chr = NA, 
      SURVIVOR.Partner1_pos = NA, 
      SURVIVOR.Partner1_gene = NA,
      SURVIVOR.Partner2_chr = NA, 
      SURVIVOR.Partner2_pos = NA, 
      SURVIVOR.Partner2_gene = NA,
      SURVIVOR.VAF = NA, 
      SURVIVOR.SV_ID = NA, 
      Found_in_SURVIVOR = NA_integer_)
    
  } else {
    
    vcf.in <- vcfR::read.vcfR(input.vcf) %>%
      vcfR::vcfR2tidy()
    
    vcf.tab <- tibble::tibble(vcf.in$fix) %>%
      tidyr::separate(ALT, into = c("ALT_base",
                                    "Partner2_chr",
                                    "Partner2_pos"),
                      remove = FALSE) %>%
      dplyr::filter(!is.na(Partner2_chr),
                    !is.na(Partner2_pos))
    
    annotated.tab <- tibble(format_SVs(vcf.tab))
    colnames(annotated.tab) <- paste("SURVIVOR", colnames(annotated.tab), sep = ".")
    
    annotated.tab <- annotated.tab %>%
      dplyr::mutate(Found_in_SURVIVOR = 1)
    
  }
  
  return(annotated.tab)
  
}

gridss_vcf_2_tibble <- function(input.vcf){
  
  vcf.in <- vcfR::read.vcfR(input.vcf) 
  
  if(nrow(vcf.in) == 0){
    
    annotated.tab <- tibble(GRIDSS.Partner1_chr = NA,
                            GRIDSS.Partner1_pos = NA, 
                            GRIDSS.Partner1_gene = NA,
                            GRIDSS.Partner2_chr = NA, 
                            GRIDSS.Partner2_pos = NA, 
                            GRIDSS.Partner2_gene = NA,
                            GRIDSS.VAF = NA, 
                            GRIDSS.SV_ID = NA, 
                            Found_in_GRIDSS = NA_integer_)
    
  } else {
    
    vcf.in <- vcfR::read.vcfR(input.vcf) %>%
      vcfR::vcfR2tidy()
    
    filename <- basename(input.vcf)
    Normal <- unlist(strsplit(filename, "--"))[2]
    
    vcf.tab <- tibble(vcf.in$gt) %>%
      dplyr::filter(!grepl(Normal, Indiv)) %>%
      dplyr::select((-c("POS", "ChromKey"))) %>%
      cbind(vcf.in$fix) %>%
      dplyr::rename(VAF = TAF) %>%
      dplyr::select(-c(starts_with("gt_"))) %>%
      tidyr::separate(ALT, 
                      sep = "[\\[\\]:]+",
                      into = c("temp1", "Partner2_chr", "Partner2_pos", "temp2"),
                      remove = FALSE) %>%
      dplyr::filter(!is.na(temp1), !is.na(Partner2_chr), !is.na(Partner2_pos), !is.na(temp2))
    
    annotated.tab <- tibble(format_SVs(vcf.tab))
    colnames(annotated.tab) <- paste("GRIDSS", colnames(annotated.tab), sep = ".")
    
    annotated.tab <- annotated.tab %>% 
      dplyr::mutate(Found_in_GRIDSS = 1)
    
  }
  return(annotated.tab)
  
  
}

manta_vcf_2_tibble <- function(input.vcf){
  
  vcf.in <- vcfR::read.vcfR(input.vcf) 
  
  if(nrow(vcf.in) == 0){
    
    annotated.tab <- tibble(Manta.Partner1_chr = NA, 
                            Manta.Partner1_pos = NA,
                            Manta.Partner1_gene = NA,
                            Manta.Partner2_chr = NA, 
                            Manta.Partner2_pos = NA, 
                            Manta.Partner2_gene = NA,
                            Manta.VAF = NA, 
                            Manta.SV_ID = NA, 
                            Found_in_Manta = NA_integer_)
    
  } else {
    
    vcf.in <- vcfR::read.vcfR(input.vcf) %>%
      vcfR::vcfR2tidy()
    
    filename <- basename(input.vcf)
    Normal <- unlist(strsplit(filename, "--"))[2]
    
    
    vcf.tab <- tibble(vcf.in$gt) %>%
      dplyr::filter(!grepl(Normal, Indiv)) %>%
      dplyr::select((-c("POS", "ChromKey"))) %>%
      cbind(vcf.in$fix) %>%
      dplyr::rename(VAF = gt_VAF) %>%
      dplyr::select(-c(starts_with("gt_"))) %>%
      tidyr::separate(ALT, into = c("ALT_base",
                                    "Partner2_chr",
                                    "Partner2_pos"),
                      remove = FALSE) %>%
      dplyr::filter(SVTYPE =="BND")
    
    annotated.tab <- tibble(format_SVs(vcf.tab))
    colnames(annotated.tab) <- paste("Manta", colnames(annotated.tab), sep = ".")
    
    annotated.tab <- annotated.tab %>%
      dplyr::mutate(Found_in_Manta = 1)
    
  }
  return(annotated.tab) #
  
}

igcaller_tsv_2_tibble <- function(input.tsv){
  
  igcaller.tab <- read.delim(input.tsv) 
  
  if(nrow(igcaller.tab) == 0){
    
    annotated.tab <- tibble(IgCaller.Partner1_chr = NA,
                            IgCaller.Partner1_pos = NA, 
                            IgCaller.Partner1_gene = NA,
                            IgCaller.Partner2_chr = NA, 
                            IgCaller.Partner2_pos = NA, 
                            IgCaller.Partner2_gene = NA, 
                            IgCaller.VAF = NA, 
                            IgCaller.SV_ID = NA, 
                            Found_in_IgCaller = NA_integer_)
    
  } else {
    
    annotated.tab <- igcaller.tab %>%
      tidyr::separate(Score, into = c("Score", "V1", "V2", "V3")) %>%
      dplyr::mutate(Score = as.numeric(Score)) %>%
      dplyr::filter(grepl("^Onco", Analysis),
                    Mechanism == "Translocation",
                    Score >= 50) %>%
      tidyr::separate(Annotation, sep = "[:; \\[\\]()]+", 
                      into = c("temp1", "temp2", "temp3", "Partner1_chr", 
                               "Partner1_pos", "Partner1_strand", "Partner2_chr", 
                               "Partner2_pos", "Partner2_strand", "Partner1_gene", 
                               "temp4", "Partner2_gene", "temp5")) %>%
      dplyr::select("Partner1_chr", "Partner1_pos", "Partner2_chr",
                    "Partner2_pos", "Partner1_gene", "Partner2_gene")
    
    colnames(annotated.tab) <- paste("IgCaller", colnames(annotated.tab), sep = ".")
    
    annotated.tab <- annotated.tab %>%
      dplyr::select(IgCaller.Partner1_chr,
                    IgCaller.Partner1_pos,
                    IgCaller.Partner1_gene,
                    IgCaller.Partner2_chr,
                    IgCaller.Partner2_pos,
                    IgCaller.Partner2_gene) %>%
      dplyr::mutate(IgCaller.SV_ID = case_when(IgCaller.Partner1_chr > IgCaller.Partner2_chr ~ substr(paste(IgCaller.Partner1_chr, IgCaller.Partner1_pos, sep = "_"), 1, 13),
                                               IgCaller.Partner2_chr > IgCaller.Partner1_chr ~ substr(paste(IgCaller.Partner2_chr, IgCaller.Partner2_pos, sep = "_"), 1, 13)),
                    Found_in_IgCaller = 1)
    
  }
  
  return(annotated.tab)
  
}


# 2. Load data
ROI.annotations <- read.csv(opt$ROI_with_IG) %>% #
  dplyr::select(-c(X)) 
ROI.annotations$names <- as.character(ROI.annotations$names)
ROI.annotations$chromosome_name <- as.character(ROI.annotations$chromosome_name)

ROI.annotations.gr <- tibble_2_GRanges(ROI.annotations$chromosome_name, 
                                       ROI.annotations$start_position, 
                                       ROI.annotations$end_position,
                                       ROI.annotations$names)


# 3. Annotate rearrangement partners
survivor.annotated.tab <- survivor_vcf_2_tibble(opt$SURVIVOR_infile) 

gridss.annotated.tab <- gridss_vcf_2_tibble(opt$GRIDSS_infile) 

manta.annotated.tab <- manta_vcf_2_tibble(opt$Manta_infile) # 

igcaller.annotated.tab <- igcaller_tsv_2_tibble(opt$IgCaller_infile) 


# 4. Merge and tidy
allSVs.merged <- survivor.annotated.tab %>%
  dplyr::full_join(gridss.annotated.tab,
                   keep = TRUE,
                   by = c("SURVIVOR.SV_ID" = "GRIDSS.SV_ID")) %>%
  dplyr::mutate(SV_ID = case_when(!is.na(SURVIVOR.SV_ID) ~ SURVIVOR.SV_ID,
                                  is.na(SURVIVOR.SV_ID) & !is.na(GRIDSS.SV_ID) ~ GRIDSS.SV_ID,
                                  is.na(SURVIVOR.SV_ID) & is.na(GRIDSS.SV_ID) ~ NA_character_)) %>%
  dplyr::full_join(manta.annotated.tab,
                   keep = TRUE,
                   by = c("SV_ID" = "Manta.SV_ID")) %>%
  dplyr::mutate(SV_ID = case_when(is.na(Manta.SV_ID) & !is.na(SV_ID) ~ SV_ID,
                                  !is.na(Manta.SV_ID) & is.na(SV_ID) ~ Manta.SV_ID,
                                  !is.na(Manta.SV_ID) & !is.na(SV_ID) ~ SV_ID,
                                  is.na(Manta.SV_ID) & is.na(SV_ID) ~ NA_character_)) %>%
  dplyr::full_join(igcaller.annotated.tab,
                   keep = TRUE,
                   by = c("SV_ID" = "IgCaller.SV_ID"))

write.csv(allSVs.merged, opt$outfile, quote = FALSE, row.names = FALSE) 


