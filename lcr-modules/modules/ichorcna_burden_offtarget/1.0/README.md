# ichorcna_burden_offtarget 1.0

## Purpose

This module derives a ploidy-free copy-number burden index from off-target
coverage. It runs the standard BAM-to-WIG preprocessing (or accepts a pre-built
WIG through `inputs.sample_wig`), runs ichorCNA once, and consumes the resulting
`correctedDepth.txt`. It does not reimplement ichorCNA normalization.

Off-target coverage does not jointly identify tumour fraction and ploidy. A
constant change in baseline ploidy can be absorbed by the model while producing
very different tumour-fraction estimates. The burden estimator therefore
median-centres corrected bin log ratios and reports their spatially correlated
copy-number structure rather than an HMM tumour-fraction point estimate. In
particular, the ploidy-3 offset `log2(3/2)` is removed exactly.

The biological limitation is important: the index measures the combination of
DNA shedding and genomic CNA complexity. CNA data alone cannot separate those
two components.

## Method

Bins are canonicalized, filtered, sorted within chromosome, and median-centred.
For non-overlapping windows of `w` bins, the estimator uses

```text
Var(window means at width w) = S + sigma^2 / w
```

A weighted linear model across `w = 1,2,4,8,16` estimates structural variance
`S` and independent technical variance `sigma^2`. `burden_amp` is
`sqrt(max(0,S))`; downstream frozen models should prefer
`burden_var_signed = S`, because zero truncation creates a floor discontinuity
for low-burden samples.

Detectability is separate from the estimate. Bins are permuted within
chromosome, preserving their marginal distribution while destroying spatial
autocorrelation. The permutation test uses the width-8 window variance, and its
null amplitude distribution supplies a per-sample detection limit. A
moving-block bootstrap provides uncertainty without destroying segment-scale
correlation.

QC flags describe assay reliability, not tumour biology. They are covariates or
sensitivity strata and are never grounds for excluding a sample.

## Inputs and execution

The module supports three input modes, in priority order:

1. Set `inputs.ichorcna_output_dir` to an existing ichorCNA `99-outputs` root
   to reuse corrected depth, CNA, segmentation, and plot files. This bypasses
   both preprocessing and ichorCNA.
2. Set `inputs.sample_wig` to an existing compiled WIG path to bypass the
   BAM-to-bigWig-to-WIG branch while still running ichorCNA.
3. Leave both values empty to generate off-target coverage from BAM and run
   ichorCNA normally.

The PON remains required for burden-output provenance. Reused ichorCNA outputs
must have been generated with the configured PON, bin size, and normalization
settings.

Run the base-R tests from this version directory:

```bash
Rscript tests/test_burden_index.R
Rscript tests/test_aggregate_burden.R
Rscript tests/test_plot_vc_diagnostic.R
```

Run the slower simulation study on demand:

```bash
Rscript tests/validate_burden_index.R
```

## Per-sample columns

All columns below are sample-intrinsic and safe to carry into a frozen
predictor, subject to pre-specifying which predictors and transformations are
used.

- `sample_id`: Sample identifier.
- `burden_index_version`: Estimator contract version.
- `burden_amp`: Zero-truncated structural amplitude, `sqrt(max(0,S))`.
- `burden_var_signed`: Untruncated structural variance `S`; preferred for models.
- `burden_amp_lo95`: Moving-block bootstrap 2.5th percentile.
- `burden_amp_hi95`: Moving-block bootstrap 97.5th percentile.
- `burden_amp_bootsd`: Bootstrap standard deviation of the amplitude.
- `detect_p`: Within-chromosome permutation p-value.
- `detect_lod95`: 95th percentile of null amplitudes.
- `detectable`: Whether `detect_p < 0.05`.
- `burden_snr`: `burden_amp / sigma_noise`.
- `sigma_noise`: Variance-components technical-noise SD.
- `sigma_noise_dlrs`: MAD of within-chromosome lag-1 differences divided by `sqrt(2)`.
- `vc_fit_r2`: Weighted variance-components fit R-squared.
- `vc_curv`: Quadratic coefficient in the weighted variance-versus-`1/w` fit; positive values are artefact-consistent curvature.
- `S_fine`: Fine-scale weighted intercept using the requested window widths.
- `S_coarse`: Coarse-scale weighted intercept using widths 4, 8, 16, and 32.
- `S_excess`: `S_fine - S_coarse`, estimating variance contributed by fine scales.
- `S_unweighted`: Fine-scale intercept with equal weights across widths.
- `amp_coarse`: `sqrt(max(0, S_coarse))`.
- `vc_V_w{w}`: Winsorised variance of window means at requested width `w`.
- `vc_nwin_w{w}`: Number of chromosome-contained windows at width `w`.
- `vc_resid_w{w}`: Residual from the weighted fine-scale fit at width `w`.
- `sd_winsor`: SD of winsorised centred bin log ratios.
- `sd_mad`: Robust MAD-based SD of centred bin log ratios.
- `fga_adaptive`: Fraction of bins beyond the adaptive noise threshold.
- `fga_fixed`: Fraction of bins with absolute centred log ratio above 0.20.
- `fga_threshold_used`: Adaptive threshold used for `fga_adaptive`.
- `n_chrom_altered`: Number of autosomes with an altered chromosome median.
- `amp_per_altered_chrom`: Burden amplitude divided by altered chromosome count.
- `burden_amp_seg`: Length-weighted SD of optional segment medians.
- `n_seg`: Number of valid segments in the optional segmentation.
- `seg_per_gb`: Valid segments per covered gigabase.
- `fga_seg`: Length fraction of segments beyond absolute centred log ratio 0.20.
- `n_bins_input`: Rows read from the bin file.
- `n_bins_used`: Bins retained after filtering.
- `n_bins_nonfinite`: Rows removed for non-finite loci or log ratios.
- `n_bins_excluded`: Bins removed by the optional exclude BED.
- `logr_offset_removed`: Median log-ratio offset removed before estimation.
- `logr_column`: Auto-detected source log-ratio column.
- `pon_id`: Basename of the normalization PON.
- `binsize_kb`: Configured bin size in kilobases.
- `winsor_p`: Symmetric winsorisation tail probability.
- `window_widths`: Comma-separated variance-component window widths.
- `detect_w`: Window width used by the detection statistic.
- `n_perm`: Number of within-chromosome permutations.
- `n_boot`: Number of moving-block bootstrap replicates.
- `block_size`: Bootstrap block size in bins.
- `seed`: Seed used only for permutation and bootstrap annotations.
- `qc_flag`: Semicolon-separated assay flags, or `pass`.

`burden_amp_seg`, `n_seg`, `seg_per_gb`, and `fga_seg` inherit dependence on
ichorCNA segmentation parameters such as `txnE` and `txnStrength`; they are
secondary descriptors.

## Cohort-relative output

Aggregation adds `burden_amp_rank`, `sigma_noise_rank`,
`cohort_sigma_noise_p50`, `cohort_sigma_noise_p90`, BH-adjusted `detect_q`, and
`detectable_fdr10`. These columns depend on the samples aggregated together and
must not be used in a frozen predictor or external validation. No data-derived
group is emitted; any categorical threshold must be fixed in the analysis plan.

The QC summary reports numeric quantiles, individual `qc_flag` counts, and the
detection rate.

## Variance-components diagnostics

The workflow automatically generates one diagnostic PDF per sample and
publishes it under
`99-outputs/{seq_type}--{genome_build}/burden_diagnostic/bin{binSize}kb/`.
The plot is part of the module `all` target.

To regenerate a diagnostic on demand from either a per-sample burden TSV or the
aggregated cohort TSV:

```bash
Rscript src/plot_vc_diagnostic.R \
  --table all_samples.burden.cohort_relative.tsv \
  --sample SAMPLE_ID \
  --out SAMPLE_ID.vc_diagnostic.pdf
```

The output extension must be `.png` or `.pdf`. The base-R script reads
`window_widths`, reconstructs the matching `vc_V_w{w}` and `vc_nwin_w{w}`
vectors, and draws the stored points against
`burden_var_signed + sigma_noise^2 / w`. It does not refit the line or modify
the burden estimate. Tables produced before the diagnostic columns were added
fail with an explicit error.

The fitted line extends to `1/w = 0` so the extrapolated structural-variance
intercept remains visible. Point size is proportional to `sqrt(vc_nwin_w{w})`,
showing why fine-window observations dominate the weighted fit. Curvature above
the line near the intercept is consistent with autocorrelated artefact inflating
`S`; curvature below the line is consistent with CNAs shorter than the wider
windows being diluted.

For cohort-level exploratory diagnostics, `analysis/vc_diagnostics.R` backfills
diagnostic columns from existing `correctedDepth.txt` files and generates
multi-sample tables and figures. Unlike the on-demand plotter, this exploratory
script uses ggplot2.

It writes sample-level and long-form TSVs, cohort overlay and residual figures,
a six-sample exemplar panel, and package/session metadata. PDF output is vector;
PNG output is 600 dpi.

```bash
Rscript analysis/vc_diagnostics.R \
  --burden all_samples.burden.cohort_relative.tsv \
  --corrected_dir 06-ichorcna_run \
  --output_dir vc_diagnostics
```

Supplying `--clinical clinical.tsv` additionally renders the mechanism figure
when the table contains `sample_id` and `PVbest_ctDNA`.
## Validation

The checked-in validation script covers amplitude recovery, noise and offset
invariance, type-I error, power and limits of detection, focal versus broad
events, and moving-block bootstrap coverage. Measured results for this
implementation are recorded here after running the full script:

```text
V1 amplitude recovery
variance_components: Pearson 0.998; Spearman 0.988; median |relative bias| 0.067
single_subtraction: Pearson 0.992; Spearman 0.937; median |relative bias| 0.109
uncorrected_sd: Pearson 0.790; Spearman 0.629; median |relative bias| 1.733

V2 noise invariance
sigma 0.05 -> recovered 0.051; burden_amp 0.076
sigma 0.08 -> recovered 0.080; burden_amp 0.077
sigma 0.12 -> recovered 0.119; burden_amp 0.077
sigma 0.15 -> recovered 0.149; burden_amp 0.078
sigma 0.20 -> recovered 0.199; burden_amp 0.078
sigma 0.25 -> recovered 0.248; burden_amp 0.080

V3 offset invariance
maximum absolute difference: 0

V4 type I error and p-value calibration
FPR 0.05: 0.067; FPR 0.10: 0.110; KS p: 0.106

V5 power and LOD
f 0.02 events 4: power 0.10; median LOD 0.019
f 0.05 events 4: power 0.05; median LOD 0.019
f 0.10 events 4: power 0.60; median LOD 0.025
f 0.20 events 4: power 1.00; median LOD 0.033
f 0.02 events 12: power 0.05; median LOD 0.018
f 0.05 events 12: power 0.45; median LOD 0.023
f 0.10 events 12: power 1.00; median LOD 0.034
f 0.20 events 12: power 1.00; median LOD 0.059
f 0.02 events 30: power 0.10; median LOD 0.019
f 0.05 events 30: power 1.00; median LOD 0.027
f 0.10 events 30: power 1.00; median LOD 0.041
f 0.20 events 30: power 1.00; median LOD 0.078

V6 focal versus broad landscapes
broad: truth 0.069; estimate 0.072
focal: truth 0.030; estimate 0.028

V7 block-bootstrap coverage
coverage: 1.000
```

## Implementation notes

- Sex chromosomes are excluded by default and included only with
  `--include_sex_chrom`.
- BED overlap uses half-open intervals (`bin_start < bed_end` and
  `bin_end > bed_start`).
- A sample with fewer than 500 retained bins receives `few_bins`; fewer than
  `min_bins` is an input error because variance-component fitting is not
  credible.
- `n_perm < 1` and `n_boot < 1` produce missing annotations rather than changing
  the deterministic burden estimate.
- Aggregation recursively discovers TSVs and intersects heterogeneous schemas
  with a warning.
