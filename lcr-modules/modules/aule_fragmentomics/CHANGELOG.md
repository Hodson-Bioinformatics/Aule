# Changelog

All notable changes to `aule_fragmentomics` are documented here.

## [Unreleased]

### Added

- DELFI-style windowed short:long profiles with sample-level LOWESS GC correction and optional healthy-reference z-scores.
- Griffin-style TFBS-centred coverage profiles, per-TF features, site-count gates, and optional healthy-reference z-scores.
- Independent healthy-control reference-table builder for DELFI bins and Griffin features.
- Regional metrics in the cohort-level sample summary.

### Fixed

- Count allele evidence from either mate and consolidate calls once per fragment.
- Require a configurable minimum base quality for allele observations.
- Declare the FASTA index as an explicit workflow dependency.

## [1.0] - 2026-08-18

### Added

- Uniform fragment filtering for deduplicated paired-end cfDNA BAMs.
- Global fragment-size distributions, short-fragment summaries, and 10 bp
  periodicity estimates.
- Reference-derived 4-mer end-motif frequencies and normalised diversity.
- Internally controlled mutant-versus-wild-type fragment-length comparisons.
- Optional application of independently trained fragment-size likelihood ratios.
- Per-sample distributions, per-variant results, and a joinable cohort table.
