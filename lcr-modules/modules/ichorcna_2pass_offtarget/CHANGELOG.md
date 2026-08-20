# Changelog

All notable changes to the `ichorcna_2pass_offtarget` module are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## Known issues (superseded by `ichorcna_burden_offtarget`)

- The Pass-2 normal grid is centred on `pass1_tf`, rather than
  `1 - pass1_tf` (`1.0/src/pass1_qc.R:246-247`); the Pass-1 TF therefore fell
  outside the Pass-2 initialisation range for 179/190 completed samples.
- `strategy` is initialized but never updated from `category`
  (`1.0/src/pass1_qc.R:282`).
- The UNDETECTABLE condition tests `t_mean_range` twice
  (`1.0/src/pass1_qc.R:286`).
- The ACCEPT gate is effectively unreachable: `loglik_gap >= 100` held for
  1/184 samples and was mutually exclusive with `t_range >= 0.20`
  (`1.0/src/pass1_qc.R:272-278`).
- The MAD registry is not wired into the Snakefile or default configuration.
- Category logic is duplicated and divergent in
  `1.0/src/final_confidence.R:258-272`.
- `scStates` and `estimateScPrevalence` are computed but ignored by the shell.
- `final_tf` is a confidence-tier-dependent mixture of Pass-1 and Pass-2 TF.
- Branches at `1.0/src/final_confidence.R:323` and in the Snakefile's
  UNRESOLVABLE arm are unreachable.
- The tests do not source the production scripts.

## [1.0] - 2026-04-28

### Added

- Two-pass ichorCNA workflow operating on off-target reads from targeted-panel BAMs.
- **Pass 1** (broad exploration): wide normal-fraction grid
  (`c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)`), all three ploidy states (2, 3, 4),
  low transition stringency (`txnE 0.99 / txnStrength 1000`) to capture
  the full likelihood landscape.
- **`pass1_qc.R`** computes five automated QC metrics from the Pass 1
  params.txt and applies decision logic:
  - `loglik_gap` — gap between top two solutions
  - `at_edge` — winner's `n_est` within 0.01 of the highest normal init
  - `t_mean_range` — spread of Student's t means (signal amplitude)
  - `convergence_fraction` — fraction of initialisations converging to same `n_est` (±5%)
  - `MAD` — GC-map correction MAD
- **Six decision outcomes**:
  | Strategy | Condition |
  |---|---|
  | `ACCEPT` | Not at edge, gap > 150, range > 0.2, MAD < 0.10, conv > 50% |
  | `UNDETECTABLE` | TF < 0.03, at edge, range < 0.08 |
  | `SMOOTH` | Noisy data (MAD > 0.10) with real signal (TF > 0.10) |
  | `RESTRICT` | Edge solution winning with multiple ploidy clusters |
  | `PLOIDY_LOCK` | Multiple clusters at similar loglik (not at edge) |
  | `UNRESOLVABLE` | None of the above |
- **Pass 2** re-runs ichorCNA with strategy-specific parameters, or forwards
  Pass 1 outputs unchanged when the strategy is `ACCEPT`, `UNDETECTABLE`,
  or `UNRESOLVABLE`.
- **`final_confidence.R`** compares both passes and assigns a final confidence
  tier: `HIGH`, `MODERATE`, `LOW`, `UNDETECTABLE`, or `UNRESOLVABLE`.
- All standard off-target BAM extraction and WIG-generation steps inherited
  from `ichorcna_offtarget 1.0`.
- Outputs include Pass 2 CNA files, Pass 1 files (for audit), per-sample
  strategy TSV, and per-sample confidence TSV.
