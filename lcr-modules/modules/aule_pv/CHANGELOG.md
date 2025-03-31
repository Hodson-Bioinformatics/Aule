# Changelog

All notable changes to the `aule_pv` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2024-11-20

This release was authored by Joanna A. Krupka.

This module performs ultra-sensitive tumour burden quantification with phased variants.

This workflow is inspired by [Kurtz et al 2021](https://www.nature.com/articles/s41587-021-00981-w). 

Input: Family collapsed bam produced by umi_collapse module.

Output: A table per patient with the total number of PV (separately for doublets and triplets), the number of informative reads (cumulative coverage for all PV rich regions at baseline), pVAF (phased variants allele frequency), estimated limit of detection.

Core functionality: 

- Input bams sorted by query name
- Find all phased variants in individual reads: iterate through paired reads, count reference mismatches per pair (ignoring indels and inserts), identify all PV combinations per read
- Sort by coordinate
- Collapse by genomic position
- Summarise per patient

