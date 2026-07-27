#! /usr/bin/env Rscript
local({
    conda <- Sys.getenv("CONDA_PREFIX", unset = NA)
    rhome <- normalizePath(R.home())
    if (!is.na(conda) && startsWith(rhome, conda)) {
        .libPaths(R.home("library"), include.site = FALSE)
    }
})
# Uses GAMBLR to remove blacklisted genes from augmented MAFs

log <- file(snakemake@log[[1]], open="wt")
sink(stdout(), type = "message")
sink(log)

library(data.table)
library(GenomicRanges)
library(rtracklayer)

genome_build        = snakemake@wildcards[["genome_build"]]
blacklist_template  = snakemake@params[["blacklist_template"]]
project_base        = snakemake@params[["project_base"]]
seq_type_param      = snakemake@params[["seq_type_blacklist"]]
drop_threshold      = snakemake@params[["threshold"]]
maf_in              = snakemake@input[["maf"]]
maf_out             = snakemake@output[["maf"]]

#### Helper function to read MAFs ####
fread_maf = function(maf_file_path){
                     maf_dt = data.table::fread(
                     file = maf_file_path,
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     verbose = FALSE,
                     data.table = TRUE,
                     showProgress = TRUE,
                     header = TRUE,
                     fill = TRUE,
                     skip = "Hugo_Symbol",
                     quote = "",
                     colClasses=c(Hugo_Symbol="character",
                                  Entrez_Gene_Id="integer",
                                  Center="character",
                                  NCBI_Build="character",
                                  Chromosome="character",
                                  Start_Position="integer",
                                  End_Position="integer",
                                  Strand="character",
                                  Variant_Classification="character",
                                  Variant_Type="character",
                                  Reference_Allele="character",
                                  Tumor_Seq_Allele1="character",
                                  Tumor_Seq_Allele2="character",
                                  dbSNP_RS="character",
                                  dbSNP_Val_Status="logical",
                                  Tumor_Sample_Barcode="character",
                                  Matched_Norm_Sample_Barcode="character",
                                  Match_Norm_Seq_Allele1="character",
                                  Match_Norm_Seq_Allele2="character",
                                  Tumor_Validation_Allele1="logical",
                                  Tumor_Validation_Allele2="logical",
                                  Match_Norm_Validation_Allele1="logical",
                                  Match_Norm_Validation_Allele2="logical",
                                  Verification_Status="logical",
                                  Validation_Status="logical",
                                  Mutation_Status="logical",
                                  Sequencing_Phase="logical",
                                  Sequence_Source="logical",
                                  Validation_Method="logical",
                                  Score="logical",
                                  BAM_File="logical",
                                  Sequencer="logical",
                                  Tumor_Sample_UUID="logical",
                                  Matched_Norm_Sample_UUID="logical",
                                  HGVSc="character",
                                  HGVSp="character",
                                  HGVSp_Short="character",
                                  Transcript_ID="character",
                                  Exon_Number="character",
                                  t_depth="integer",
                                  t_ref_count="integer",
                                  t_alt_count="integer",
                                  n_depth="integer",
                                  n_ref_count="integer",
                                  n_alt_count="integer",
                                  all_effects="character",
                                  Allele="character",
                                  Gene="character",
                                  Feature="character",
                                  Feature_type="character",
                                  Consequence="character",
                                  cDNA_position="character",
                                  CDS_position="character",
                                  Protein_position="character",
                                  Amino_acids="character",
                                  Codons="character",
                                  Existing_variation="character",
                                  ALLELE_NUM="integer",
                                  DISTANCE="logical",
                                  STRAND_VEP="integer",
                                  SYMBOL="character",
                                  SYMBOL_SOURCE="character",
                                  HGNC_ID="character",
                                  BIOTYPE="character",
                                  CANONICAL="character",
                                  CCDS="character",
                                  ENSP="character",
                                  SWISSPROT="character",
                                  TREMBL="character",
                                  UNIPARC="character",
                                  RefSeq="character",
                                  SIFT="character",
                                  PolyPhen="character",
                                  EXON="character",
                                  INTRON="character",
                                  DOMAINS="character",
                                  AF="character",
                                  AFR_AF="character",
                                  AMR_AF="character",
                                  ASN_AF="logical",
                                  EAS_AF="character",
                                  EUR_AF="character",
                                  SAS_AF="character",
                                  AA_AF="character",
                                  EA_AF="character",
                                  CLIN_SIG="character",
                                  SOMATIC="character",
                                  PUBMED="numeric",
                                  MOTIF_NAME="logical",
                                  MOTIF_POS="logical",
                                  HIGH_INF_POS="logical",
                                  MOTIF_SCORE_CHANGE="logical",
                                  IMPACT="character",
                                  PICK="integer",
                                  VARIANT_CLASS="character",
                                  TSL="integer",
                                  HGVS_OFFSET="integer",
                                  PHENO="character",
                                  MINIMISED="logical",
                                  GENE_PHENO="integer",
                                  FILTER="character",
                                  flanking_bps="character",
                                  vcf_id="character",
                                  vcf_qual="character",
                                  gnomAD_AF="logical",
                                  gnomAD_AFR_AF="logical",
                                  gnomAD_AMR_AF="logical",
                                  gnomAD_ASJ_AF="logical",
                                  gnomAD_EAS_AF="logical",
                                  gnomAD_FIN_AF="logical",
                                  gnomAD_NFE_AF="logical",
                                  gnomAD_OTH_AF="logical",
                                  gnomAD_SAS_AF="logical",
                                  vcf_pos="integer",
                                  gnomADg_AF="character"))

  return(maf_dt)
}

seq_type = if (seq_type_param == "both") c("genome", "capture") else snakemake@wildcards[["seq_type"]]

this_maf = fread_maf(maf_in)
blacklist <- rtracklayer::import(blacklist_template)  # GRanges

if (nrow(this_maf)>0){
   ## ---- build GRanges from MAF ----
    gr_maf <- GRanges(
    seqnames = this_maf$Chromosome,
    ranges   = IRanges(
        start = this_maf$Start_Position,
        end   = this_maf$End_Position
        )
    )

    ## ---- overlap test ----
    this_maf$blacklisted <- as.integer(overlapsAny(gr_maf, blacklist))

    ## ---- write output ----
    fwrite(
        this_maf,
        file  = maf_out,
        sep   = "\t",
        quote = FALSE,
        na    = ""
    )

} else {
    message(paste("WARNING: Detected 0 vartiants for the sample", snakemake@wildcards[["tumour_id"]]))
    fwrite(
        this_maf,
        file  = maf_out,
        sep   = "\t",
        quote = FALSE,
        na    = ""
    )
}

sink()
