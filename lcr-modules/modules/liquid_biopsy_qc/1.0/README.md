# liquid_biopsy_qc 1.0

Assay QC for liquid-biopsy ichorCNA and burden-index analyses. The module
records metrics and advisory flags; it never filters samples and never uses
detectability as a QC criterion.

## Assay and coverage modes

Set `options.assay_mode` to one of:

- `capture`: remove every alignment overlapping `options.panel_bed`, then run
  samtools stats/flagstat and mosdepth on the temporary off-target BAM.
- `swgs`: run samtools stats/flagstat and mosdepth directly on the full BAM.
  Capture-only Picard HS metrics and off-target aliases are recorded as `NA`
  with `not_applicable_swgs`.

With the default `compute_coverage: True`, coverage is generated inside this
module. Set it to `False` to reuse four persistent, precomputed sidecars via
`inputs.offtarget_stats`, `offtarget_flagstat`,
`offtarget_mosdepth_summary`, and `offtarget_mosdepth_regions`. Those legacy
input names are retained for compatibility; in sWGS mode they represent
genome-wide rather than off-target coverage.

Coverage is summarized in generic columns:

- `coverage_scope`: `off_target` or `genome_wide`
- `coverage_reads` and `coverage_read_frac_observed`
- `mean_coverage_depth`, `median_coverage_depth`, and
  `coverage_depth_cv`

Capture runs additionally populate the backward-compatible
`offtarget_*` columns and compare observed off-target fraction with Picard HS
metrics. sWGS runs do not apply that comparison.

The mosdepth window is controlled by `coverage_bin_size_bp`; the default is
1 Mb. `min_mean_coverage_depth` is assay-specific and must be set from the
validated capture or sWGS design rather than estimated from the analysed cohort.

## Minimal configuration

Capture:

```yaml
options:
    assay_mode: capture
    compute_coverage: True
    panel_bed: data/panel_extended.bed
    coverage_bin_size_bp: 1000000
```

sWGS:

```yaml
options:
    assay_mode: swgs
    compute_coverage: True
    panel_bed: ""
    coverage_bin_size_bp: 1000000
    thresholds:
        min_mean_coverage_depth: <validated-sWGS-threshold>
```

The ichorCNA corrected depth, parameters, PON, GC and mappability inputs must
come from the same assay design. In particular, a capture off-target PON must
not be applied to sWGS data.

Set `skip_bam_metrics: True` only to collect ichorCNA/burden tiers without
full-BAM or coverage metrics. Missing values are written as `NA` and explained
in `missing_reason`.

## SLURM grouping

Coverage, BAM metrics, and per-sample collection belong to the connected
`liquid_biopsy_qc_sample` group. For each sample they execute sequentially in
one SLURM allocation:

```text
coverage -> BAM metrics -> collection
```

This produces one longer scheduler submission per sample rather than three
shorter submissions. Because the rules are sequential, Snakemake requests the
maximum CPU and memory needed by any stage instead of summing their resources.
No `--group-components` option is required. Independent samples remain
independent group jobs and are scheduled normally with `--jobs`.
