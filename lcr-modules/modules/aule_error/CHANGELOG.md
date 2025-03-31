# Changelog

All notable changes to the `aule_error` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2024-08-02

This release was authored by Ilias Moutsopoulos.

This module can be used to estimate error rates for each panel position, either from healthy samples to be used as background, or from tumour samples to see the differences. 

The errors considered are: UMI errors, where a UMI “family” is discordant between its members, R1/R2 error, where the R1 and R2 reads are overlapping and discordant at a location, and directional errors, where the base at a position does not match the reference. 

The rates for all error types are calculated by a  python script using pysam: A bam file containing UMI collapsed reads produced by fgbio CallMolecularConsensusReads is used as input, aligned and sorted by query name. For each read pair, the aligned genomic positions are extracted, as well as the per-position consensus depth and consensus errors (cd & ce tags respectively), and the sequence agreement for the positions that the R1 and R2 reads are overlapping. The above metrics are then tallied for every position in the panel as the reads in the bam file are parsed.

Finally, the metrics are aggregated across all samples to produce a final table, by panel position. When using normal samples, this is conceptually equivalent to an "error panel of normals" that can be used as a baseline for error filtering. Tables for individual samples are also reported.
