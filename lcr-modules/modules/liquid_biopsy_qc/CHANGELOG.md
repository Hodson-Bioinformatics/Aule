# Changelog

## 1.0

- Initial liquid-biopsy QC implementation.
- Add native capture off-target coverage generation with bedtools, samtools, and mosdepth.
- Add explicit capture/sWGS assay modes and generic coverage metrics.
- Retain external coverage-sidecar mode for compatibility with ichorcna_offtarget.
- Run coverage, BAM metrics, and QC collection sequentially in one SLURM group
  per sample, with Picard heap headroom and measured per-stage resources.
