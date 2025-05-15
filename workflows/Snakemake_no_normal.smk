#!/usr/bin/env snakemake

'''
This Snakefile is made to run all the modules compatible with Capture workflow.
Compatibility of a workflow can be checked by referring to the pairing_config parameter 
present in a default.yaml file of that module.
'''
##### SETUP #####

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")
aule_github_clone = "/home/jak75/src/Aule" # Path to installed AULE 

##### REFERENCE_FILES WORKFLOW #####

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../lcr-modules/workflows/reference_files/2.4/reference_files.smk" 
    configfile:
        aule_github_clone +  "/lcr-modules/workflows/reference_files/2.4/config/default.yaml"

##### CONFIGURATION FILES #####

# Load module-specific configuration
configfile: aule_github_clone + "/lcr-modules/modules/aule_preprocessing/1.0/config/default.yaml"
configfile: aule_github_clone + "/lcr-modules/modules/mutect2/2.0/config/default.yaml"
configfile: aule_github_clone + "/lcr-modules/modules/vcf2maf/1.3/config/default.yaml"

# Load project-specific config, which includes the shared
# configuration and some module-specific config updates

configfile: "config_demo.yaml"
config["lcr-modules"]["_shared"]["samples"] = SAMPLES

##### MODULE SNAKEFILES #####

# Load module-specific snakefiles
include: aule_github_clone +"/lcr-modules/modules/aule_preprocessing/1.0/aule_preprocessing.smk"
include: aule_github_clone +"/lcr-modules/modules/mutect2/2.0/mutect2.smk"
include: aule_github_clone + "/lcr-modules/modules/vcf2maf/1.3/vcf2maf.smk"

##### TARGETS ######

rule aule_test:
    input:
        rules._aule_preprocessing_all.input,
        rules._mutect2_all.input,
        rules._vcf2maf_all.input
