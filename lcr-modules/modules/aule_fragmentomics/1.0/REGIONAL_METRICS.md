# DELFI and Griffin regional metrics

The regional collector implements measured profiles and reference application as
separate operations. It never trains on the analysed tumour cohort and does not
label a profile as a cancer probability.

## Configuration

Enable either analysis independently:

```yaml
lcr-modules:
    aule_fragmentomics:
        inputs:
            # Download/provide this raw BED4/5 catalogue. The workflow cleans it.
            griffin_raw_sites_bed: "data/dlbcl_tfbs.raw.bed.gz"
        options:
            assay_partition: "off_target"
            enable_delfi: True
            enable_griffin: True
            generate_delfi_bins: True
            prepare_griffin_sites: True
```

With `generate_delfi_bins: True`, the module creates autosomal, non-overlapping
5 Mb windows from the FASTA index, calculates GC once from the FASTA, annotates
blacklist overlap, and removes windows above `delfi_max_blacklist_fraction`.
The BED6 is retained under `99-outputs/assets/delfi_bins.bed`. Set
`generate_delfi_bins: False` to supply `delfi_bins_bed` instead; BED column 4
may contain a stable bin ID and column 5 precomputed GC.

`griffin_sites_bed` is BED4/5: chromosome, start, end, TF/site-set name, and an
optional site ID. Site coordinates must use the same genome build and contig
style as the BAM and FASTA. At least `griffin_min_sites_per_tf` sites (default
100) are required before a feature can receive `pass` status or a reference
z-score; thousands of sites per TF are preferable at low depth.

With `prepare_griffin_sites: True`, `griffin_raw_sites_bed` is checked against
the FASTA index, edge- and blacklist-filtered, deduplicated, optionally restricted
by `griffin_tf_allowlist`, and sorted. The BED5 is retained under
`99-outputs/assets/griffin_sites.bed`. Set `prepare_griffin_sites: False` to
use an already prepared `griffin_sites_bed`.

The reference TSVs are optional. Without them, raw and GC-corrected profiles are
still produced, but normalized scores remain `NA` and the summary states that a
reference was not supplied. References must be made from independent healthy
controls matched for collection, processing, library chemistry, assay partition,
panel, and probe lot where applicable.

## DELFI calculation

Each proper-pair fragment is counted once from its positive TLEN alignment after
the module's uniform MAPQ, flag, length, and blacklist filters. Its midpoint is
assigned to one DELFI bin. Per bin the collector reports short (100–150 bp) and
long (151–220 bp) counts and:

```text
raw_log2_short_long = log2((short + pseudocount) / (long + pseudocount))
```

A local-linear LOWESS fit estimates the sample-level relationship between bin GC
fraction and the raw log ratio. The corrected value preserves the sample median:

```text
corrected = raw - LOWESS(GC) + median(raw)
```

When a healthy reference is supplied, each corrected bin is transformed using
its independently estimated `reference_mean` and `reference_sd`. The module emits
per-bin z-scores, mean absolute z, and RMS z. These are fragmentation-deviation
summaries, not the proprietary/trained DELFI cancer classifier.

DELFI remains partly collinear with copy number. Downstream modelling should test
incremental value over ichorCNA tumour fraction/CNA burden or regress matched
bin-level copy-number measurements before constructing a classifier.

## Griffin-style calculation

For every TF/site set, fragment coverage is accumulated in configurable bins
across the default ±1 kb site-centred window. Coverage is divided by the number of
sites so TFs with different catalogue sizes remain comparable. Reference GC is
calculated for each profile bin. A local-linear LOWESS model of coverage against
mean site GC is fitted per sample and TF, followed by multiplicative correction
and within-profile mean normalization.

The per-TF features are:

- `central_coverage`: mean normalized coverage within ±30 bp.
- `mean_coverage`: mean raw fragments covering a profile bin per site.
- `amplitude`: max–min normalized coverage within ±250 bp.

This is an auditable Griffin-style site-centred implementation, not a wrapper
around the external Griffin software package. The exact GC model and feature
definitions are recorded in the output columns and configuration. Healthy
reference z-scores are applied only to TF/feature combinations with a positive
reference SD and a passing site-count/GC-correction gate.

## Reference schemas and construction

DELFI reference:

```text
bin_id  reference_mean  reference_sd  n_reference_samples
```

Griffin reference:

```text
tf  feature  reference_mean  reference_sd  n_reference_samples
```

Build both from completed healthy-control outputs using:

```bash
python src/build_regional_reference.py \
  --delfi-inputs control1.delfi.tsv,control2.delfi.tsv \
  --griffin-inputs control1.griffin_features.tsv,control2.griffin_features.tsv \
  --delfi-output delfi_healthy_reference.tsv \
  --griffin-output griffin_healthy_reference.tsv \
  --min-samples 10
```

The builder uses sample-level values, calculates the sample standard deviation,
drops invariant features, and refuses duplicate sample/key observations.

The workflow can run the builder automatically, but it never guesses which
samples are controls:

```yaml
        inputs:
            delfi_reference_inputs: "/controls/{sample_id}.delfi.tsv"
            griffin_reference_inputs: "/controls/{sample_id}.griffin_features.tsv"
        options:
            build_delfi_reference: True
            build_griffin_reference: True
            reference_sample_ids: [CONTROL01, CONTROL02, CONTROL03]
            reference_min_samples: 3
```

Production references should normally contain at least 10 matched healthy
controls. Generated references are retained in `99-outputs/assets/`.

## What remains external

The module generates deterministic derivatives locally. You still provide the
reference FASTA and FAI, a genome-build-matched blacklist, BAM/BAI files, and the
biological TF-site catalogue used by Griffin. Variant VCFs are needed only for
allele-specific metrics. A healthy reference is optional; when requested, its
completed control result files and sample IDs must be supplied explicitly.

## Outputs

- `delfi/*.delfi.tsv`: counts, GC model, corrected ratios, and per-bin z-scores.
- `griffin_profiles/*.griffin_profile.tsv`: complete TF-centred coverage profiles.
- `griffin_features/*.griffin_features.tsv`: three features and optional z-scores
  per TF.
- `per_sample/*.regional_fragmentomics.tsv`: compact DELFI/Griffin status and
  deviation summaries.
- `aggregate/cohort_fragmentomics.tsv`: global, allele, DELFI, and Griffin scalar
  summaries joined by sample ID.
