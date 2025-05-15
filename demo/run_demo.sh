#!/bin/bash

nice snakemake \
  --use-conda \
  --conda-frontend conda \
  -s Snakemake_demo.smk \
  --jobs 1 \
  aule_test

#### Consider adding these lines for HPC tests
 # --cluster-cancel scancel \
 #  --cluster 'sbatch --mem={resources.mem_mb} --ntasks {threads} --partition={resources.slurm_partition} --time={resources.runtime}' \
 #  --default-resources slurm_partition="icelake" slurm_account="hodson-sl3-cpu" runtime=60 \
 #  --slurm \