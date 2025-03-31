# Changelog

All notable changes to the `manta` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2020-07-17

This release was authored by Bruno Grande.

- Update the sample IDs in the VCF header before parsing the VCF file with cyvcf2 in case the BAM files don't have sample IDs, causing the `somaticSV` to have two `N/A` sample columns.

## [1.0]

This release was authored by Bruno Grande.
