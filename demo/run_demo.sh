#!/bin/bash

nice snakemake \
  --use-conda \
  -s Snakemake_demo.smk \
  --dry-run \
  aule_test