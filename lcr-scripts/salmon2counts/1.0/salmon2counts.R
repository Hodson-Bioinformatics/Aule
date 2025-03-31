### Salmon to Counts Matrix ### 
# Uses tximport to convert salmon quant.sf files 
# to transcript- and gene-level counts matrices. 

#!/usr/bin/env Rscript
#
# Usage: 
#   Rscript salmon2counts.R <path/to/salmon/results> <path/to/transcripts.gtf> <path/to/output> <sample_table.tsv> [<sample_id_column>]
#
# Notes: 
#   This script is intended for use with the Salmon-1.0 module in LCR-modules. 
#   It expects to find the salmon quant.sf files at the input path, following 
#   the pattern {sample_id}.quants.sf. These files should be in the 
#   99-outputs subdirectory of the salmon-1.0 output directory. 
#
#   The transcripts.gtf (or .gff) file should be the same set of transcripts used to 
#   run Salmon. 
#
#   The sample table should adhere to LCR-modules guidelines. If no sample_id
#   column is specified, the script will look for the values to fill in for 
#   {sample_id}.quant.sf in the `sample_id` column. The sample table should at 
#   least be subset to samples with seq_type == "mrna", but can be further subset
#   to samples of interest prior to running this script. Only samples in the 
#   samples table will be included in the output tables. 

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(tximport)
  library(GenomicFeatures)
})

# Determine arguments -----------------------------------------------------

# Parse command-line arguments (when running non-interactively)
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("salmon_path", "transcript_gtf", "output_path", "sample_table", "sample_id_column")
args <- setNames(args, arg_names[1:length(args)])

# Args for testing interactively
# args <- list()
# args$salmon_path <- "/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/gambl/salmon-1.0/99-outputs/quant_to_hg38/mrna/"
# args$sample_table <- "/projects/rmorin/projects/gambl-repos/gambl-lhilton/data/metadata/gambl_samples_available.tsv"
# args$transcript_gtf <- "/projects/rmorin_scratch/lcr-modules-references-DEV/genomes/hg38/annotations/gencode_annotation-33.gtf"
# args$output_path <- "~/salmon_test"

# Print args for debugging
print(paste("salmon_path:",args$salmon_path,"sample_table:",args$sample_table,"transcript_gtf:",args$transcript_gtf,"output_path",args$output_path))

if(is.null(args$sample_id_column)) args$sample_id_column <- "sample_id"

# Load gtf file and convert to txdb or use existing one from output directory-------------------
tx2gene_file = paste0(args$output_path,"/tx2gene.tsv")
if(file.exists(tx2gene_file)){
  message("Loading existing tx2gene file")
  tx2gene<-read.table(tx2gene_file)
}else{
  message("Loading transcript file and generating tx2gene...")
  txdb <- makeTxDbFromGFF(args$transcript_gtf)
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- select(txdb, k, "GENEID", "TXNAME")
  fwrite(tx2gene, tx2gene_file, sep = "\t")
}

# Load sample table and get sample IDs and Salmon paths-------------------

message("Loading sample table and finding available data...")
metadata <- read_tsv(args$sample_table) %>% filter(seq_type == "mrna")
sample_table <- metadata[args$sample_id_column]
colnames(sample_table) <- "sample_id"

sample_table <- sample_table %>% 
  mutate(salmon_path = str_c(args$salmon_path, "/", sample_id, ".quant.sf"), 
         salmon_available = file.exists(salmon_path)) 

salmon_missing <- sample_table[!sample_table$salmon_available, "sample_id", drop = TRUE]
if(length(salmon_missing) != 0){
  warning("The following samples do not have Salmon results: ", 
        paste0(salmon_missing, collapse = ", "))
}

  
sample_table <- sample_table[sample_table$salmon_available, ]  

# Import Salmon files with tximport --------------------------------------

salmon_files <- sample_table$salmon_path
names(salmon_files) <- sample_table$sample_id

message("Loading ", length(salmon_files), " samples with tximport...")

# Transcript-level counts
transcripts <- list()
transcripts$txi <-
  tximport(salmon_files,
           type = "salmon",
           txIn = TRUE,
           txOut = TRUE)

transcripts$counts <-
  round(transcripts$txi$counts) %>%
  as.data.frame() %>%
  rownames_to_column("transcript")

transcripts$lengths <-
  round(transcripts$txi$length) %>%
  as.data.frame() %>%
  rownames_to_column("transcript")

# Summarize transcripts to genes
genes <- list()
genes$txi <- summarizeToGene(transcripts$txi, tx2gene = tx2gene)

genes$counts <-
  round(genes$txi$counts) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

genes$lengths <-
  round(genes$txi$length) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

# Output --------------------------------------------------------------------------------------

message("Writing transcript-level counts to file...")
fwrite(transcripts$counts, paste0(args$output_path, "/salmon.transcripts.counts.tsv"), sep = "\t")
fwrite(transcripts$lengths, paste0(args$output_path, "/salmon.transcripts.lengths.tsv"), sep = "\t")
write_rds(transcripts$txi, paste0(args$output_path, "/salmon.transcripts.txi.rds"))

message("Writing gene-level counts to file...")
fwrite(genes$counts, paste0(args$output_path, "/salmon.genes.counts.tsv"), sep = "\t")
fwrite(genes$lengths, paste0(args$output_path, "/salmon.genes.lengths.tsv"), sep = "\t")
write_rds(genes$txi, paste0(args$output_path, "/salmon.genes.txi.rds"))

metadata_out <- metadata[metadata$sample_id %in% sample_table$sample_id, ] 

write_tsv(metadata_out, paste0(args$output_path, "/salmon.metadata.tsv"))

message("Complete.")

