# Changelog

All notable changes to the `ichorcna_offtarget` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-11-17

### Changed

- Default `options.bamCoverage` is now
  `--ignoreDuplicates --extendReads --samFlagExclude 2048 --minMappingQuality 20`
  (was `--ignoreDuplicates --extendReads`), matching the Aule-setups
  `ichorcna-offtarget` setup. Supplementary alignments and reads with MAPQ < 20
  are no longer counted. Runs that relied on the old default will get slightly
  different coverage; set `options.bamCoverage` in your config to keep the old
  behaviour. Coverage should be counted with the same filters as the PoN: the
  12-control Novogene PoN was built with `--samFlagExclude 1028`.

## [1.0] - 2026-04-23

This release was authored by Joanna A. Krupka.

This module is a variant of ichorcna/2.0 that first excludes on-target panel regions from the input BAMs using `bedtools intersect -v` before computing coverage. 
This is useful for off-target copy number analysis from targeted sequencing (e.g., DIRECT panel).
A panel regions BED file (typically extended by 2 kb and collapsed) must be supplied as `panel_bed` in inputs.
The off-target BAM is created as a temporary file and is not retained after the run.

### Added

- Optional `inputs.sample_wig` support. A non-empty wildcard path bypasses the
  BAM-to-WIG branch; an empty value retains the existing preprocessing behavior.
- Automatic cohort aggregation of per-sample ichorCNA parameter files into
  `99-outputs/aggregate/all_samples.ichorcna_summary.tsv`.
- Persistent samtools, flagstat, and 1 Mb mosdepth sidecars are emitted while
  the temporary off-target BAM exists for downstream liquid-biopsy QC.
