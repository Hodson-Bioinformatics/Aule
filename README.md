# Aule

![Lymphoma-3](https://github.com/Hodson-Bioinformatics/Aule/assets/32344189/b063b1cc-ca76-454d-b803-ac39c8818a94)

## What is AULE?
Aule is the Automated Ultrasensitive Lymphoma Evaluation pipeline performing end-to-end analysis of circulating tumour DNA (ctDNA) in patients with Large B-cell Lymphoma (LBCL). 
AULE has a modular architecture following the framework of the [LCR pipeline](https://lcr-modules.readthedocs.io/en/latest/), 
developed by Prof. Ryan Morin's lab. 

## Quick start 

1. Download a copy of the repository to your local directory 

```
git clone https://github.com/Hodson-Bioinformatics/Aule.git
```

2. Create a baseline conda environment from `setup/envs/aule_py3.7.yaml`

```
cd Aule
conda env create -f setup/envs/aule_py3.7.yaml
```
