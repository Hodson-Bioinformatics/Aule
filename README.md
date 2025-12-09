# AULE

![aule_logo_v1](https://github.com/user-attachments/assets/f53e6b12-0ede-4310-84cf-2633fa547592)

**Aule** (Automated Ultrasensitive Lymphoma Evaluation) is a pipeline performing end-to-end analysis of circulating tumour DNA (ctDNA) in patients with Large B-cell Lymphoma (LBCL). 
AULE has a modular architecture following the framework of the [LCR pipeline](https://lcr-modules.readthedocs.io/en/latest/), developed by Prof. Ryan Morin's lab. AULE can be used as a comprehensive standalone software package or integrated as an add-on to an existing LCR pipeline.
 
## 🔬 Background

Recent studies have generated considerable interest in the potential of ctDNA to transform cancer treatment and facilitate precision medicine trials with treatments tailored to individual risk profiles and biological subtypes. However, the majority of current ctDNA analysis platforms are proprietary "black-box" assays that offer limited opportunities for customisation, further development, and rigorous validation. To address these limitations, we have developed a comprehensive open-source pipeline that renders ctDNA analysis both accessible and cost-effective for research.

## ⭐️ Why AULE? 

AULE addresses critical challenges in ctDNA analysis for lymphoma monitoring:

- **End-to-end workflow**  - From raw FASTQs to Minimal Residual Disease (MRD) statistics in one modular pipeline
- **Ultra-sensitivity**  - Detect single nucleotide variants (SNVs) at frequencies as low as 0.01% using UMIs-based error modelling
- **Tumour-specific tracking**  - Monitor MRD status with phased variants unique to each patient (as pioneered by [Kurtz et al. 2021](https://www.nature.com/articles/s41587-021-00981-w))
- **Open science**  - Fully reproducible with documented methods, version control and open codebase.
- **Modular design** - Flexible execution of individual modules or end-to-end workflows, adaptable to already existing workflows

## 📒 AULE usage 

Learn more on [AULE wiki](https://github.com/Hodson-Bioinformatics/Aule/wiki) 

#### ⚡️ Quick links

[Installation](https://github.com/Hodson-Bioinformatics/Aule/wiki/Installation)  
[AULE architecture and an example workflow](https://github.com/Hodson-Bioinformatics/Aule/wiki/Architecture-and-basic-workflow)  
[How to setup your first AULE run](https://github.com/Hodson-Bioinformatics/Aule/wiki/Typical-run-setup)  

## 🎯 AULE in action

[DIRECT trial](https://cctu.org.uk/portfolio/cancer/trials-closed-to-recruitment-in-follow-up/direct) - a multisite prospective molecular profiling trial of patients undergoing first-line therapy for aggressive B-cell Non-Hodgkin lymphoma, see our [preprint](https://doi.org/10.1101/2025.04.14.25325806) on MedRxiv

## 💪 Contributors 

[Joanna A. Krupka](https://github.com/ashakru)  
[Ilias Moutsopoulos](https://github.com/emouts)  
[Natasha H Cutmore](https://github.com/nuttynutmore)   
[Daniel J Hodson](https://www.stemcells.cam.ac.uk/people/pi/hodson)  
