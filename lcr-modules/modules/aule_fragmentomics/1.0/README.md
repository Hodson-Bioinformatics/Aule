# aule_fragmentomics 1.0

cfDNA fragmentomics for deduplicated paired-end BAMs, with optional independent healthy-reference application.
The module implements the highest-priority metrics that can be computed without
inventing a reference cohort or region catalogue:

1. Global fragment-size distribution and summary metrics.
2. Reference-derived 4-mer fragment-end motifs and normalised Shannon diversity.
3. Mutant-versus-wild-type fragment-length comparisons at supplied biallelic SNVs.
4. Optional application of independently trained fragment-size likelihood ratios.
5. DELFI-style windowed short:long profiles with GC correction and optional healthy-reference z-scores.
6. Griffin-style TFBS-centred nucleosome profiles with per-TF features and optional healthy-reference z-scores.

The module describes samples and does not exclude them or train a classifier.

## Uniform fragment definition

A fragment is counted from the alignment with positive `TLEN`, which counts each
pair once without read-name bookkeeping. It must be a proper pair, pass
`MAPQ >= 30`, have none of SAM flags `3852`, have length 50–250 bp, and not
overlap the configured ENCODE blacklist. All cut-offs are explicit config values.

The two 5-prime end motifs are the reference 4-mer at the left fragment boundary
and the reverse complement of the reference 4-mer ending at the right boundary.
All 256 motifs, including zero-count motifs, are emitted.

## Allele-specific analysis

The variant VCF analysis currently accepts PASS biallelic SNVs. Query-to-reference
coordinates are resolved through the CIGAR alignment. Results are emitted for
every variant, including named statuses for unsupported or underpowered records.
Either mate may supply the allele observation; mates are consolidated by read
name so each fragment is counted once. Allele bases must pass the configured
minimum base quality (30 by default).
The sample-level distribution counts each fragment name once; fragments with
conflicting calls across variant sites are reported and omitted from that summary.

Size weighting is application-only. `size_weights` must contain `length` and
`likelihood_ratio` (or `weight`) columns and must have been trained on independent
high-tumour-fraction samples. Missing length weights produce `NA` weighted results
with `incomplete_weight_table`; they are never silently replaced.

## Outputs

- `99-outputs/per_sample/*.fragmentomics.tsv`: joinable scalar summaries.
- `99-outputs/fragment_size/*.fragment_size_hist.tsv`: complete size histogram.
- `99-outputs/end_motifs/*.end_motifs.tsv`: 256 motif counts and frequencies.
- `99-outputs/allele_fragment_size/*.allele_fragment_size_hist.tsv`: mutant and
  wild-type size distributions.
- `99-outputs/per_variant/*.per_variant.tsv`: counts, shifts, and two-sample KS
  results for each variant.
- `99-outputs/delfi/*.delfi.tsv`: windowed DELFI counts, corrected ratios, and z-scores.
- `99-outputs/griffin_profiles/*.griffin_profile.tsv`: TF-centred coverage profiles.
- `99-outputs/griffin_features/*.griffin_features.tsv`: per-TF coverage features.
- `99-outputs/assets/delfi_bins.bed`: locally generated, GC-annotated DELFI windows.
- `99-outputs/assets/griffin_sites.bed`: validated and filtered Griffin site catalogue.
- `99-outputs/assets/*_healthy_reference.tsv`: optional references built only
  from explicitly listed healthy-control outputs.
- `99-outputs/aggregate/cohort_fragmentomics.tsv`: combined sample table.

Run the module separately on full, on-target, or off-target BAMs when those assay
partitions are required. `options.assay_partition` is required and is copied into
every scalar summary so capture-derived and sWGS distributions are not presented
as directly interchangeable.

## Regional metrics

See [REGIONAL_METRICS.md](REGIONAL_METRICS.md) for input schemas, exact DELFI and Griffin-style calculations, reference construction, validation gates, and interpretation limits.

## Deliberately deferred

OCF, WPS/phasing, promoter metrics, preferred-end discovery, and
FDI are not included in version 1.0. They require one or more assets not supplied
with this module: matched healthy reference cohorts, GC models, TFBS/accessibility
sets, probe-boundary annotations, promoter designs, or validated profile tensors.
Adding placeholder implementations would create plausible but uninterpretable
features. Jagged-end inference is excluded because the required library chemistry
is absent.
