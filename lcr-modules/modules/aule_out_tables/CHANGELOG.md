# Changelog

All notable changes to the `aule_out_tables` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2] - 2025-01-28

This release was authored by Joanna Krupka

The module was reworked to focus on mutational data. An optional filtering step was also added, to filter plasma mutations using the annotations created in aule_plasma.

## [1.1] - 2024-08-29

This release was authored by Natasha Cutmore.

Alteration to structure of metadata to account for genes having separate OG/TSG and SHM status.
Removed structural variant data - will summarise in a separate module.
Tidied nomenclature for hotspot mutations to match output from VEP105 using canonical TXIDs (except for TCF3 and ETS1).

## [1.0] - 2023-02-27

This release was authored by Ilias Moutsopoulos.

This module is responsible for creating combined tables by patient, by sample and by mutation. The module is designed to be run at the end of the pipeline and bring together all data for a study.

Mutational data in particular is annotated with a lymphoma whitelist for commonly occuring mutations, and drivers are annotated through a set of rules curated through the available literature.
