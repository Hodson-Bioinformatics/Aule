# Changelog

All notable changes to the `ichorcna_pon_offtarget` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-11-17

### Changed

- Default `options.bamCoverage.opts` now uses `--samFlagExclude 2048` (was
  `1028`): `--ignoreDuplicates --extendReads --samFlagExclude 2048 --minMappingQuality 20`.
  Supplementary alignments are no longer counted; unmapped reads are still
  dropped by `--minMappingQuality 20` and duplicates by `--ignoreDuplicates`.
  This matches the `ichorcna_offtarget` default and the Aule-setups convention,
  so a PoN built with defaults pairs with samples run with defaults. PoNs built
  before this change (e.g. the 12-control Novogene PoN) used `1028`; set
  `options.bamCoverage.opts` to rebuild one of them exactly.

## [1.0] - 2026-04-23

This release was authored by Joanna A. Krupka.

This module is a variant of ichorcna_pon/1.0 that first excludes on-target panel regions from the input BAMs using `bedtools intersect -v` before computing coverage for panel-of-normals creation.
This is useful for building off-target PONs from targeted sequencing data (e.g., DIRECT panel).
A panel regions BED file (typically extended by 2 kb and collapsed) must be supplied as `panel_bed` in inputs.
The off-target BAM is created as a temporary file and is not retained after the run.
