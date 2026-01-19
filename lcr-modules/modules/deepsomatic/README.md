# DeepSomatic Module

## Overview

The DeepSomatic module integrates Google's DeepSomatic variant caller into the LCR-modules framework. DeepSomatic is a deep learning-based somatic variant caller that can operate in both matched tumor-normal and tumor-only modes.

## Features

- **Matched tumor-normal mode**: Standard somatic variant calling with paired tumor and normal samples
- **Tumor-only mode**: Somatic variant calling using only tumor samples with germline filtering via population databases
- **Docker integration**: Uses official Google DeepSomatic Docker containers
- **Multiple sequencing types**: Support for WGS, WES, and RNA-seq data
- **FFPE optimization**: Specialized models for FFPE (Formalin-Fixed Paraffin-Embedded) samples
- **Configurable PON filtering**: Built-in Panel of Normals filtering for tumor-only mode
- **Multi-threading support**: Parallel processing for improved performance

## Supported Model Types

### Matched Tumor-Normal Mode
- `FFPE_WGS`: Whole genome sequencing of FFPE samples
- `FFPE_WES`: Whole exome sequencing of FFPE samples

### Tumor-Only Mode  
- `FFPE_WGS_TUMOR_ONLY`: Whole genome sequencing tumor-only mode
- `FFPE_WES_TUMOR_ONLY`: Whole exome sequencing tumor-only mode

## Configuration

### Required Inputs

For matched tumor-normal mode:
```yaml
inputs:
  sample_bam: "path/to/tumor.bam"
  sample_bai: "path/to/tumor.bam.bai"
  normal_bam: "path/to/normal.bam"
  normal_bai: "path/to/normal.bam.bai"
```

For tumor-only mode:
```yaml
inputs:
  sample_bam: "path/to/tumor.bam"
  sample_bai: "path/to/tumor.bam.bai"
  normal_bam: ""
  normal_bai: ""
```

### Options

```yaml
options:
  # Docker image version
  docker_image: "google/deepsomatic:1.9.0"
  
  # Enable default PON filtering for tumor-only mode
  use_default_pon_filtering: "true"
  
  # Optional: restrict analysis to specific regions
  regions:
    capture: "chr1,chr2,chr3"  # Example
    
  # Additional command-line arguments
  deepsomatic_matched: "--min_mapping_quality=20"
  deepsomatic_tumor_only: "--min_base_quality=10"
```

## Resource Requirements

The module is configured for high-performance computing environments:

- **Memory**: 32 GB RAM (configurable)
- **CPU**: 8 threads (configurable)
- **Disk**: 50 GB temporary storage (configurable)
- **Docker**: Requires Docker installation and permissions

## Output Files

The module generates the following outputs:

1. **VCF files**: Compressed variant call files (`.vcf.gz`)
2. **Index files**: Tabix index files (`.vcf.gz.tbi`)
3. **Logs**: Detailed execution logs in the logs directory
4. **Intermediate results**: Optional intermediate files for debugging

Output naming convention:
- Matched mode: `{sample_id}.matched.deepsomatic.vcf.gz`
- Tumor-only mode: `{sample_id}.tumor_only.deepsomatic.vcf.gz`

## Usage Examples

### Basic Matched Tumor-Normal Analysis

```yaml
lcr-modules:
  deepsomatic:
    inputs:
      sample_bam: "results/bam/{sample_id}.tumor.bam"
      sample_bai: "results/bam/{sample_id}.tumor.bam.bai"
      normal_bam: "results/bam/{sample_id}.normal.bam"
      normal_bai: "results/bam/{sample_id}.normal.bam.bai"
```

### Tumor-Only Analysis with Custom Regions

```yaml
lcr-modules:
  deepsomatic:
    inputs:
      sample_bam: "results/bam/{sample_id}.tumor.bam"
      sample_bai: "results/bam/{sample_id}.tumor.bam.bai"
      normal_bam: ""
      normal_bai: ""
    options:
      regions:
        capture: "chr1:1000000-2000000"
      use_default_pon_filtering: "true"
```

## Dependencies

- Docker (for running DeepSomatic containers)
- Snakemake
- Reference genome files (FASTA + index)

## References

- [DeepSomatic GitHub Repository](https://github.com/google/deepsomatic)
- [DeepSomatic Documentation](https://github.com/google/deepsomatic/tree/r1.9/docs)
- [Original DeepSomatic Paper](https://doi.org/10.1101/2022.06.02.494516)