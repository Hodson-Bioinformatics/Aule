# Changelog

All notable changes to `ichorcna_burden_offtarget` are documented here.

## [1.0] - 2026-07-28

### Added

- One-pass off-target ichorCNA preprocessing with a ploidy-free burden index.
- Per-sample variance-components burden, detectability, uncertainty, complexity,
  and assay-QC annotations.
- Cohort aggregation with explicitly cohort-relative ranks, quantiles, FDR, and
  QC summaries.
- Base-R production tests and an on-demand simulation validation suite.
- Dynamically named `vc_V_w{w}` and `vc_nwin_w{w}` columns preserving the
  variance-component observations for default and custom window widths.
- A base-R, single-sample variance-components diagnostic plot generated solely
  from stored burden-table values.
- Automatic per-sample diagnostic PDF generation and publication through the
  Snakemake `all` target.
- Optional `inputs.ichorcna_output_dir` support to bypass preprocessing and
  ichorCNA by reusing compatible published corrected-depth and segmentation
  outputs.
