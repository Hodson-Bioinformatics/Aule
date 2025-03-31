# Changelog

All notable changes to the `aule_plasma` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2025-01-06

This release was authored by Ilias Moutsopoulos.

The module was reworked to start from the 2+ variant caller output vcf of the slms3 module. It annotates variants with various metrics that can be used for filtering, as well as an extra flag of whether a variant appears in the 3+ vcf (and is therefore called with higher confidence). Variants can then be filtered using these metrics, either manually or through aule_out_tables.

## [1.0] - 2024-05-16

This release was authored by Joanna A. Krupka.

The module takes as input the union vcf file from slms3 and annotates each variant with a variety of metrics that can be used for filtering. As this is effectively a de novo variant caller and would require extensive testing to ensure its robustness, it was not finalised in 1.0 but remained in development.
