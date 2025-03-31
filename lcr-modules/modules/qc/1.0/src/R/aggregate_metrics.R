#!/usr/bin/env Rscript
#

suppressWarnings(suppressPackageStartupMessages({

  message("Loading packages...")

  required_packages <- c(
    "data.table",
    "tidyverse"
  )

  for (pkg in required_packages){
    library(pkg, character.only = TRUE)
  }
}))



# Collect the base quality
base_scores = fread(
  snakemake@input[[2]],
  skip = 7
) %>%
  as.data.frame

QualityScore =
  base_scores %>%
  mutate(
    totalQ = as.numeric(QUALITY) * as.numeric(COUNT_OF_Q),
    totalQ = sum(totalQ),
    totalN = sum(COUNT_OF_Q)
  ) %>%
  mutate(AverageBaseQuality = round(totalQ / totalN, 2)) %>%
  distinct(AverageBaseQuality) %>%
  `rownames<-`(snakemake@wildcards$sample_id)


# handle samtools metrics
command = paste("grep ^SN", snakemake@input[[1]], "| cut -f 2-")
samtools_scores = fread(cmd = command,fill=TRUE,sep = "\t") %>%
  select(-V3) %>%
  mutate(V1=gsub(":", "", V1)) %>%
  column_to_rownames("V1") %>%
  `names<-`(snakemake@wildcards$sample_id) %>%
  t %>%
  as.data.frame

SamtoolsMetrics =
  samtools_scores %>%
  select(
    `insert size average`,
    `average length`,
    `pairs on different chromosomes`,
    `raw total sequences`,
    `reads mapped`,
    `reads unmapped`,
    `reads duplicated`
  )

SamtoolsMetrics =
  SamtoolsMetrics %>%
  mutate(
    ProportionReadsDuplicated = round(`reads duplicated` / `raw total sequences`, 2),
    ProportionReadsMapped = round(`reads mapped` / `raw total sequences`, 2)
  )

colnames(SamtoolsMetrics) = c(
  "AverageInsertSize",
  "AverageReadLength",
  "PairsOnDiffCHR",
  "TotalReads",
  "TotalUniquelyMapped",
  "TotalUnmappedreads",
  "TotalDuplicatedreads",
  "ProportionReadsDuplicated",
  "ProportionReadsMapped"
)

# handle coverage metrics
if (snakemake@wildcards$seq_type == "capture"){
  command = paste("grep BAIT_SET", snakemake@input[[3]] ,"-A 1")
} else {
  command = paste("grep GENOME_TERRITORY", snakemake@input[[3]] , "-A 1")
}

CoverageMetrics = fread(cmd=command) %>%
  as.data.frame %>%
  `rownames<-`(snakemake@wildcards$sample_id)

if (snakemake@wildcards$seq_type == "capture"){
  CoverageMetrics =
    CoverageMetrics %>%
    select(MEAN_TARGET_COVERAGE,
           ZERO_CVG_TARGETS_PCT,
           PCT_TARGET_BASES_10X,
           PCT_TARGET_BASES_30X,
           PCT_TARGET_BASES_50X,
           PCT_TARGET_BASES_100X,
           PCT_TARGET_BASES_250X,
           PCT_TARGET_BASES_500X,
           PCT_TARGET_BASES_1000X,
           PCT_TARGET_BASES_2500X,
           PCT_TARGET_BASES_5000X,
           PCT_TARGET_BASES_10000X,
           PCT_TARGET_BASES_100000X)
} else {
  CoverageMetrics =
    CoverageMetrics %>%
    select(MEAN_COVERAGE,
           PCT_10X,
           PCT_30X,
           PCT_50x,
           PCT_100x,
           PCT_250x,
           PCT_500x,
           PCT_1000x,
           PCT_2500x,
           PCT_5000x,
           PCT_10000x,
           PCT_100000x) %>%
    mutate(ProportionTargetsNoCoverage=NA, .after = MEAN_COVERAGE)
}

colnames(CoverageMetrics) = c(
  "MeanCorrectedCoverage",
  "ProportionTargetsNoCoverage",
  "ProportionCoverage10x",
  "ProportionCoverage30x",
  "ProportionCoverage50x",
  "ProportionCoverage100x",
  "ProportionCoverage250x",
  "ProportionCoverage500x",
  "ProportionCoverage1000x",
  "ProportionCoverage2500x",
  "ProportionCoverage5000x",
  "ProportionCoverage10000x",
  "ProportionCoverage100000x"
)

outputMetrics = cbind(SamtoolsMetrics, CoverageMetrics, QualityScore) %>%
  relocate(AverageBaseQuality)

outputMetrics$SeqType = snakemake@wildcards$seq_type

outputMetrics =
  outputMetrics %>%
  relocate(SeqType) %>%
  as.data.frame %>%
  rownames_to_column("UID")

write_tsv(outputMetrics, snakemake@output[[1]])
