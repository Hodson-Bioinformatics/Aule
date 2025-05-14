#!/usr/bin/env snakemake

'''
This Snakefile is made to run all the modules compatible with Capture workflow.
Compatibility of a workflow can be checked by referring to the pairing_config parameter 
present in a default.yaml file of that module.
'''
##### SETUP #####

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

##### REFERENCE_FILES WORKFLOW #####

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../lcr-modules/workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../lcr-modules/workflows/reference_files/2.4/config/default.yaml"

##### CONFIGURATION FILES #####

# Load module-specific configuration
configfile: "../lcr-modules/modules/aule_preprocessing/1.0/config/default.yaml"

# Load project-specific config, which includes the shared
# configuration and some module-specific config updates

configfile: "config_demo.yaml"
config["lcr-modules"]["_shared"]["samples"] = SAMPLES

##### MODULE SNAKEFILES #####

# Load module-specific snakefiles
include: "../lcr-modules/modules/aule_preprocessing/1.0/aule_preprocessing.smk"

##### TARGETS ######

rule aule_test:
    input:
        rules._aule_preprocessing_all.input
