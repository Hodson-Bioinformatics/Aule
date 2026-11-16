# Changelog

All notable changes to `ichorcna_burden_qc` are documented here. Before
2026-11-16 the module was called `ichorcna_burden_offtarget`.

## [1.0] - 2026-11-16 (renamed from `ichorcna_burden_offtarget` 1.0)

### Changed

- Renamed `ichorcna_burden_offtarget` to `ichorcna_burden_qc`: module directory,
  Snakefile, config key (`lcr-modules: ichorcna_burden_qc`) and every rule
  (`_ichorcna_burden_qc_*`). Output directories are now
  `ichorcna_burden_qc-1.0/` with subdirectories `00-inputs`, `01-burden`,
  `02-burden_plot` and `99-outputs`.
- `inputs.ichorcna_output_dir` is now required: the module reads an existing
  `ichorcna_offtarget` fit and stops at startup if it is unset.
- The R scripts run in the conda env `r`, which points at the same
  `ichorcna_offtarget` environment as before (`conda_envs.ichorcna_run`).
- Default `mem_mb` for `burden` and `aggregate` lowered from 4000 to 3000, so a
  one-core job fits a ~3.4 GB-per-core Slurm partition.
- `99-outputs/` now holds only this module's results (`burden/`,
  `burden_diagnostic/`, `aggregate/`). The ichorCNA files (corrected depth,
  CNA, segments, plot) stay in `ichorcna_offtarget`'s `99-outputs/` instead of
  being linked here again.

### Removed

- Computing off-target coverage from BAMs and running ichorCNA inside this
  module (`inputs.sample_bam`, `sample_bai`, `panel_bed`, `sample_wig`;
  `options.bamCoverage`, `ichorcna_branch`, `ichorcna_run`; the rules
  `install_ichorcna`, `input_bam`, `exclude_panel`, `bamCoverage`,
  `bigwigToWig`, `convert_wig`, `compile_wigs` and `ichorcna_run`). Their
  defaults differed from `ichorcna_offtarget` (e.g. no `--samFlagExclude` or
  `--minMappingQuality` in `bamCoverage`, `--ploidy 'c(2,3,4)'`), so a fit made
  here could silently disagree with the PON and with the fit QC checks. The
  install rule also clashed with `ichorcna_offtarget`'s, forcing a `ruleorder`
  in any Snakefile that included both.

### Unchanged

- `src/burden_index.R`, `src/aggregate_burden.R`, `src/plot_vc_diagnostic.R`
  and all `options.burden` settings: given the same fit, the burden table is
  the same.

## [1.0] - 2026-07-28 (as `ichorcna_burden_offtarget`)

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
