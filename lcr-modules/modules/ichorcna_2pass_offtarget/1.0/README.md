# ichorCNA 2-Pass Workflow — R Scripts Architecture

Generated: 2026-05-06 (Updated to match Snakemake workflow)

## Overview

This directory contains the R scripts that implement the automated 2-pass ichorCNA analysis workflow. The Snakemake module (`ichorcna_2pass_offtarget.smk`) orchestrates the full pipeline via conda environments; these scripts execute the statistical decision logic and QC assessment at each stage.

## Directory Structure

```
ichorcna_2pass_offtarget/1.0/
├── ichorcna_2pass_offtarget.smk    # Snakemake workflow definition
├── src/                             # R scripts (this directory)
│   ├── pass1_qc.R                  # Script 1: Pass 1 QC metrics & strategy
│   ├── final_confidence.R          # Script 2: Per-sample confidence assignment
│   ├── aggregate_params.R          # Script 3: Cohort-level params merge
│   ├── mad_registry.R              # Script 4: MAD threshold management
│   ├── combine_plots.R             # Script 5: PDF combination
│   └── combine_genome_wide_plots.R # Legacy (alternative implementation)
├── tests/                           # Test suite
│   ├── test_pass1_qc.R            # Tests for pass1_qc.R
│   ├── test_final_confidence.R    # Tests for final_confidence.R
│   ├── test_mad_registry.R        # Tests for mad_registry.R
│   └── test_conv_quality.R        # Convergence quality tests
└── README.md                        # This file
```

## Scripts Overview

### Script 1: `pass1_qc.R` — Pass 1 QC and Strategy Decision

**Purpose:** Execute per-sample after Pass 1, compute QC metrics, and assign Pass 2 strategy.

**Invoked by:** Rule `_ichorcna_2pass_offtarget_pass1_qc`

**CLI Interface:**
```bash
Rscript pass1_qc.R \
  --params <Pass 1 .params.txt> \
  --strategy <output strategy.tsv> \
  --sample <sample_id> \
  --max_normal_init <0.99> \
  [--mad_registry_file <path>] \
  [--mad_registry_run_type off_target|all_target]
```

**Input:** 
- ichorCNA Pass 1 `.params.txt` file (header block + solutions table)
- Optional: MAD registry CSV for threshold lookup

**QC Metrics Computed:**
- `p1_MAD`: GC-Map correction from params file
- `p1_t_mean_range`: Max − min of Student's t mean values
- `p1_loglik_gap`: Difference between top 2 solutions by log-likelihood
- `p1_convergence`: Fraction of solutions with n_est within ±0.05 of winner
- `p1_at_edge`: Boolean, winner n_est close to max_normal_init (≥0.98)
- `p1_multiple_ploidy`: Boolean, two ploidy clusters with gap < 50 in loglik
- `p1_ploidy_loglik_gap`: Gap between top two ploidy clusters

**Pass 1 Accept Check:**
```
p1_accepted = (
  !p1_at_edge AND
  p1_loglik_gap ≥ 100 AND
  p1_t_mean_range ≥ 0.20 AND
  p1_MAD < mad_p90_threshold AND
  p1_convergence ≥ 0.60
)
```

**Strategy Categories (Priority Order):**
1. **ACCEPT** — p1_accepted override (all other conditions ignored)
2. **UNDETECTABLE** — t_range < 0.08 AND pass1_tf < 0.03
3. **PLOIDY_AMBIGUOUS** — p1_multiple_ploidy AND gap < 50
4. **NOISY_WEAK** — t_range < 0.15 AND MAD > mad_p90
5. **WEAK_SIGNAL** — t_range < 0.15 AND MAD ≤ mad_p90
6. **NOISY_SIGNAL** — t_range ≥ 0.15 AND MAD > mad_p90
7. **CLEAN_SIGNAL** — default (t_range ≥ 0.15 AND MAD ≤ mad_p90)

**Pass 2 Parameters by Category:**

| Category | txnE | txnStrength | maxCN | Normal Grid | Ploidy | scStates | estimateScPrevalence |
|----------|------|-------------|-------|-------------|--------|----------|----------------------|
| ACCEPT | — | — | — | — | — | — | — |
| UNDETECTABLE | 0.9999 | 10000 | 3 | P1±15% | c(2) | c() | FALSE |
| NOISY_WEAK | 0.9999 | 10000 | 3 | P1±15% | c(2) | c() | FALSE |
| WEAK_SIGNAL | 0.999 | 3000 | 4 | P1±20% | c(P1) | c() | FALSE |
| CLEAN_SIGNAL | 0.999 | 3000 | 5 | P1±15% | c(P1) | c() | FALSE |
| NOISY_SIGNAL | **0.999** | 5000 | 5 | P1±20% | c(P1) | c() | FALSE |
| PLOIDY_AMBIGUOUS | 0.99 | 1000 | 7 | broad | c(dominant) | c(1,3) | TRUE |

**Critical Note:** NOISY_SIGNAL uses `txnE = 0.999` (NOT 0.9999), which preserves model distinction on real-signal samples.

**Output:** `strategy.tsv` (tab-separated key-value pairs, no header)

**MAD Threshold Resolution:**
- If `--mad_registry_file` provided AND registry has ≥10 samples for run_type: load P90 from registry
- Otherwise: use hardcoded fallback (off_target=0.10, all_target=0.13) + warning

---

### Script 2: `final_confidence.R` — Per-Sample Confidence Tier Assignment

**Purpose:** Read Pass 1 & 2 params, strategy TSV; assign final confidence tier and estimate flag.

**Invoked by:** Rule `_ichorcna_2pass_offtarget_confidence`

**CLI Interface:**
```bash
Rscript final_confidence.R \
  --pass1_params <Pass 1 .params.txt> \
  --pass2_params <Pass 2 .params.txt> \
  --strategy <strategy.tsv> \
  --output <confidence.tsv> \
  --sample <sample_id>
```

**Input:**
- Pass 1 and Pass 2 `.params.txt` files
- Strategy TSV from `pass1_qc.R`

**Pass 2 Metrics Computed:**
- Parsed identically to Pass 1 (loglik_gap, convergence, at_edge, multiple_ploidy)
- `tf_abs_diff = |pass1_tf − pass2_tf|`
- `conv_quality`: convergence quality label (see below)

**Convergence Quality Assessment (Applied in order, first match wins):**

| Label | Convergence | Loglik Gap |
|-------|-------------|-----------|
| TRIVIAL | ≥0.8 | <5 |
| HOLLOW | ≥0.6 | <15 |
| STRONG | ≥0.6 | ≥30 |
| GAP_DRIVEN | ≥0.3 | ≥80 |
| GENUINE | ≥0.4 | ≥15 |
| POOR | — | — |

**Special Case: ACCEPT/UNDETECTABLE Strategies**
When Pass 2 outputs are copies of Pass 1:
- Set all p2_* metrics equal to p1_* metrics
- Set tf_abs_diff = 0
- Set conv_quality = "STRONG"

**Confidence Tier Routing (Applied in order):**

1. **UNDETECTABLE strategy** → HIGH, "undetectable_confirmed"
2. **ACCEPT strategy** → HIGH, "pass1_accepted"
3. **Pass 2 failed entirely** (p2_at_edge AND p2_convergence ≤0.3) → UNRESOLVABLE
4. **HIGH via Pass 2:**
   - Category-specific gap/convergence thresholds
   - !p2_at_edge AND tf_abs_diff ≤0.05 AND conv_quality ∈ {STRONG, GAP_DRIVEN}
   - gap ≥ category_threshold AND convergence ≥ category_threshold
   - Flag: "pass1_pass2_agree"
5. **MODERATE:**
   - !p2_at_edge AND p2_convergence ≥0.3
   - Flag: "pass2_stable" (if tf_agree) or "pass2_diverged_from_pass1"
6. **LOW (default):**
   - If p2_convergence ≥0.3 AND !p2_at_edge: "unresolvable_pass2_best_effort"
   - Elif p2_at_edge: "pass2_at_edge"
   - Else: "pass2_low_convergence"

**Post-Assignment Suffix Flags:**

- **INDETERMINATE_HIGH:** Applied if tier==LOW AND final_tf≥0.10 AND p2_loglik_gap<10 AND p2_convergence<0.40
  - Suffix appended to flag
- **CONVERGED_FLAT:** Applied if conv_quality==TRIVIAL AND strategy ∈ {CLEAN_SIGNAL, NOISY_SIGNAL}
  - Suffix appended to flag
  - If tier==HIGH: demote to MODERATE (parameter artefact)

**Output:** `confidence.tsv` (tab-separated, one row per sample, with header)

---

### Script 3: `aggregate_params.R` — Merge All Samples' Solutions Tables

**Purpose:** Combine per-sample Pass 1 and Pass 2 solutions into a single cohort-level TSV.

**Invoked by:** Rule `_ichorcna_2pass_offtarget_aggregate_params`

**CLI Interface:**
```bash
Rscript aggregate_params.R <output.tsv> <p1_params...> <p2_params...>
```

**Positional Arguments:**
- `arg[1]`: output TSV path
- `arg[2..N+1]`: Pass 1 `.params.txt` files (N samples)
- `arg[N+2..2N+1]`: Pass 2 `.params.txt` files (same sample order)

**Processing:**
- For each sample, extract solutions table from both passes
- Add columns:
  - `sample_id`: extracted from filename
  - `pass`: 1 or 2
  - `is_winner`: TRUE for highest loglik per sample+pass
- Merge all rows and write to single TSV

**Output:** `all_samples_all_models.tsv` — one row per solution per sample per pass

---

### Script 4: `mad_registry.R` — MAD Threshold Management

**Purpose:** Manage MAD audit trail across batches; compute P90 and Tukey fence thresholds.

**Function Signature:**
```r
update_mad_registry(
  confidence_tsv,
  outdir,
  mad_registry_path = NULL,
  update_registry = FALSE,
  batch_id = NULL,
  run_types = c("off_target", "all_target")
)
```

**Three Operational Modes:**

**Mode A — Compute from current batch (default, no registry):**
- Compute P90 and Tukey fence from confidence_tsv$p1_MAD
- Write new registry and thresholds
- Log: "MAD thresholds computed from N samples (Mode A)"

**Mode B — Use supplied thresholds (registry provided, update=FALSE):**
- Load existing registry, extract thresholds
- Append current batch rows with threshold_source="user_supplied"
- Do NOT recompute thresholds
- Log: "Using pre-computed MAD thresholds from <path> (Mode B). Thresholds not recomputed."

**Mode C — Merge and recompute (registry provided, update=TRUE):**
- Load existing registry, append current batch
- Recompute thresholds from combined set
- Log before/after sample counts per run_type
- Log: "Merged and recomputed (Mode C)"

**Output Files (written to `outdir/`):**

1. **`ichorCNA_mad_registry.csv`** — persistent audit trail (append-only)
   - Columns: sample_id, run_type, p1_MAD, batch_id, date_added, threshold_source
   - One row per sample per batch

2. **`ichorCNA_mad_thresholds.csv`** — threshold summary
   - Columns: run_type, n_samples, mad_median, mad_p75, mad_p90, mad_tukey, threshold_source, registry_path, computed_at
   - One row per run_type

**Integration with pass1_qc.R:**
- `pass1_qc.R` looks for existing registry to load thresholds per run_type
- If registry has <10 samples: fall back to hardcoded thresholds + warning
- If thresholds change materially (Δ > 0.005), rerunning per-sample QC may change categories

---

### Script 5: `combine_plots.R` — Combine Pass 1 and Pass 2 PDFs

**Purpose:** Merge Pass 1 and Pass 2 genome-wide PDFs into a single document for QC review.

**Invoked by:** Rule `_ichorcna_2pass_offtarget_aggregate_pdf`

**CLI Interface:**
```bash
Rscript combine_plots.R <output.pdf> <p1_plots...> <p2_plots...>
```

**Positional Arguments:**
- `arg[1]`: output PDF path
- `arg[2..N+1]`: Pass 1 PDFs (N samples)
- `arg[N+2..2N+1]`: Pass 2 PDFs (same sample order)

**Processing:**
- Validate all PDF files exist
- Combine using `pdftools::pdf_combine()` (alternating P1 then P2)
- Fallback to concatenation if combining not available

**Output:** Multi-page PDF with alternating Pass 1 and Pass 2 plots

---

## Workflow Sequence

```
1. BAM → WIG conversion (Snakemake: bamCoverage, bigWigToWig, etc.)
           ↓
2. Pass 1 ichorCNA run (Snakemake: runIchorCNA.R with broad initial grid)
           ↓
3. pass1_qc.R → strategy.tsv (QC check, category routing, P2 parameters)
           ↓
4. Pass 2 ichorCNA run (Snakemake: conditional on strategy)
           ↓
5. final_confidence.R → confidence.tsv (tier assignment, estimate flags)
           ↓
6. aggregate_params.R → all_samples_all_models.tsv (merge across samples)
           ↓
7. mad_registry.R → MAD thresholds (update registry, log threshold changes)
           ↓
8. combine_plots.R → combined PDF (QC gallery)
```

## Testing

Run the testthat suite:
```bash
Rscript -e "library(testthat); test_dir('tests/')"
```

Or per-script:
```bash
Rscript -e "library(testthat); source('tests/test_pass1_qc.R')"
```

## Key Constraints

1. **Format Strictness:** Strategy TSV must use format `key\tvalue` (no spaces), as Snakemake reads via `grep -P "^key\t" | cut -f2`.

2. **NOISY_SIGNAL txnE:** Must be 0.999, NOT 0.9999. The latter collapses likelihood on real-signal samples (trivial convergence artefact).

3. **CONVERGED_FLAT Demotion:** Samples with TRIVIAL convergence quality (e.g., parameter overfitting) in CLEAN_SIGNAL/NOISY_SIGNAL categories must be demoted from HIGH to MODERATE confidence.

4. **MAD Thresholds:** Never hardcoded in final output except as labelled fallbacks. Always log which source was used (registry vs. fallback).

5. **Confidence TSV Structure:** Exactly one header row + one data row per sample. Aggregate rule concatenates assuming this.

6. **File Paths:** All scripts use `file.path()` for path construction, never `paste()`.

## Dependencies

- **R packages:** optparse, pdftools, testthat
- **Conda environments:** ichorcna_run, aggregate_plots, etc. (managed by Snakemake)

## Future Extensions

- Automated rerunning of per-sample QC if MAD thresholds change materially
- Batch-level quality dashboards (currently per-sample only)
- Integration with external structural variant callers (DELLY, GRIDSS, etc.)
