# Changelog

All notable changes to the `deepsomatic` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-01-06

This release was authored by Joanna A. Krupka.

### Added
- Initial implementation of DeepSomatic variant calling module
- Support for both matched tumor-normal and tumor-only modes
- Integration with Docker containers for DeepSomatic v1.9.0
- Configurable model types for different sequencing platforms (WGS, WES)
- Support for FFPE samples with specialized model types
- Optional PON (Panel of Normals) filtering for tumor-only mode
- Configurable regions for targeted analysis
- Multi-threading support for improved performance
- Comprehensive resource and memory management