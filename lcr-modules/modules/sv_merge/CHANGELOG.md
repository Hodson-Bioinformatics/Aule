# Changelog

All notable changes to the `direct_out_tables` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2024-12-24

This release was authored by Natasha Cutmore.

Significant updates to module:
Removed use of IgCaller output.
Removed use of SURVIVOR.

Restructured module to take output from SV callers and 1) merge all SVs 2) filter this file for SVs only in regions of interest (ROI) 3) filter this on read support and summarise into per-sample table output.

Intersection of SV caller output now uses StructuralVariantAnnotation package from R (rather than SURVIVOR).

## [1.0] - 2023-05-30

This release was authored by Natasha Cutmore and Ilias Moutsopoulos.

The purpose of this module is to merge the output of multiple structural variant callers, conceptually similarly to how slms3 merges the output of multiple somatic variant callers. The module takes the output of gridss, manta, and igcaller and merges them into a single table. The SURVIVOR package is used to merge the vcfs of gridss and manta, and a custom R script is used to merge the resulting file with the igcaller output (as it is a tsv and not a vcf). Contrary to slms3, variants are not intersected or filtered based on their presence in multiple callers. The module outputs a single table with all variants, annotated with the caller(s) that identified them.
