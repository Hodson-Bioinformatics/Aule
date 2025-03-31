# Changelog

All notable changes to the `aule_umi_collapse` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0] - 2024-06-26

This release was authored by Joanna A. Krupka.

Added CallOverlappingConsensusBases to identify R1/R2 discordant bases. Also added some intermediate sorting steps as the some commands in fgbio are very particular about receiving sorted input.


## [1.0] - 2024-05-16

This release was authored by Joanna A. Krupka.

This module collapses UMIs using the fgbio worflow (GroupReadsByUmi followed by CallMolecularConsensusReads). The GATK workflow used in aule_preprocessing preserves all reads while storing UMI information in the BAM file, which can be incompatible with further analysis steps. The BAM file is realigned after UMI collapsing to ensure that the reads (collapsed UMIs) are correctly aligned to the reference genome. A queryname-sorted BAM is returned for compatibility with the phased variant module.
