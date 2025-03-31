#!/usr/bin/env Rscript

# Load packages and functions -----------------------------------------------------------------

library(biomaRt)
library(magrittr)
library(data.table)
library(tidyverse)

`%?%` <- function(x, y) {
  if (length(x) == 0) {
    return(y)
  } else if (length(x) == 1) {
    if (is.na(x) | is.null(x)) {
      return(y)
    } else {
      return(x)
    }
  } else {
    ifelse(is.na(x) | is.null(x), y, x)
  }
}

# Argument parsing ----------------------------------------------------------------------------

args      <- commandArgs(trailingOnly = TRUE)
maf_path  <- args[1] %?% "results/strelka/8-tabulate_canonical/BLGSP.strelka.dna_aug.rna_aug.qss_filt.canonical.maf"
out_path  <- args[2] %?% "results/strelka/9-tabulate_noncanonical/BLGSP.better_noncanonical_transcripts.tsv"
min_diff  <- as.integer(args[3]) %?% 2L
ensv      <- as.integer(args[4]) %?% 86L

vep_cols <- c("gene", "consequence", "change", "transcript", "refseq")

# Load and tidy data --------------------------------------------------------------------------

maf <- fread(maf_path)


# Extract mutations -------------------------------------------------------

mutations <- maf[
  !is.na(all_effects),
  .(Tumor_Sample_Barcode, Transcript_ID, all_effects)
  ][, id := .I] %>%
  setkey(id)

mutated_tx <-
  mutations[, .(id, split = strsplit(all_effects, ";", fixed = TRUE))] %>%
  unnest() %>%
  separate(split, vep_cols, ",", extra = "merge") %>%
  as.data.table() %>%
  setkey(id)

nonsyn_cons <- c(
  "coding_sequence_variant", "conservative_missense_variant",
  "disruptive_inframe_deletion", "disruptive_inframe_insertion",
  "exon_loss_variant", "frameshift_variant", "inframe_deletion",
  "inframe_insertion", "initiator_codon_variant", "missense_variant",
  "protein_altering_variant", "rare_amino_acid_variant",
  "splice_acceptor_variant", "splice_donor_variant", "start_lost",
  "stop_gained", "stop_lost", "transcript_ablation")

mutated_tx <- merge(mutated_tx, mutations)[,all_effects := NULL]
mutated_tx[, `:=`(selected = transcript == Transcript_ID,
                  nonsyn = consequence %in% nonsyn_cons)]
mutated_tx <- mutated_tx[, .(sample = Tumor_Sample_Barcode, transcript,
                             gene, nonsyn, selected, refseq)]
mutated_tx <- unique(mutated_tx) %>% setkey("transcript")


# Get transcript lengths ----------------------------------------------------------------------

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                      version = ensv)
if (ensv >= 78) {
  attrs <- c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length",
             "external_gene_source")
  filters <- list(transcript_biotype = "protein_coding")
  transcripts <-
    getBM(attributes = attrs, filters = names(filters), values = filters,
          mart = ensembl, uniqueRows = TRUE) %>%
    filter(!grepl("Clone-based", external_gene_source)) %>%
    select(gene_id = ensembl_gene_id, transcript = ensembl_transcript_id,
           length = transcript_length) %>%
    as.data.table() %>%
    setkey("transcript")
} else {
  attrs <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
             "exon_chrom_start", "exon_chrom_end")
  transcripts <-
    getBM(attributes = attrs, mart = ensembl, uniqueRows = TRUE) %>%
    rename(gene_id = ensembl_gene_id, transcript = ensembl_transcript_id) %>%
    mutate(exon_length = exon_chrom_end - exon_chrom_start) %>%
    group_by(gene_id, transcript) %>%
    summarise(length = sum(exon_length)) %>%
    as.data.table() %>%
    setkey("transcript")
}


 # Incorporate transcript length ---------------------------------------------------------------

prioritized <- merge(mutated_tx, transcripts)
prioritized <- prioritized[, num_nonsyn := length(sample[nonsyn]),
                           by = .(gene_id, gene, transcript, length, selected, refseq)]
prioritized <- prioritized[, num_nonsyn_selected := as.integer(max(0, num_nonsyn[selected])),
                           by = .(gene_id)]
prioritized <- prioritized[!selected & refseq != ""]
prioritized <- prioritized[, -c("sample", "nonsyn", "selected")] %>% unique()
prioritized <- prioritized[, diff := as.integer(num_nonsyn - num_nonsyn_selected),
                           by = .(gene_id)]
prioritized <- prioritized[order(gene_id, -diff, -length)]
prioritized <- prioritized[, .SD[1], by = .(gene_id)]
prioritized <- prioritized[diff >= min_diff]
prioritized <- prioritized[order(-diff),
                           .(transcript, refseq, gene_id, gene, num_nonsyn_selected,
                             num_nonsyn, diff)]

readr::write_tsv(prioritized, out_path)
