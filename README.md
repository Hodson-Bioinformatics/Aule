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

### 📒 AULE usage 

Learn more on [AULE wiki](https://github.com/Hodson-Bioinformatics/Aule/wiki) 

### 🎯 AULE in action

AULE has been used in the following research studies:

- **Krupka et al. (2025)** - [*The DIRECT study: A roadmap for ctDNA-based risk prediction, molecular profiling and MRD detection in Diffuse Large B Cell Lymphoma*](https://www.medrxiv.org/content/10.1101/2025.04.14.25325806v1) - The primary publication, where AULE was developed and applied to a large cohort of LBCL patients in the DIRECT clinical trial 

### 📖 Citing AULE:

If you use AULE in your research, please cite:

```
   @article{krupka2025direct, 
   year     = {2025}, 
   title    = {The {DIRECT} study: A roadmap for {ctDNA}-based risk prediction, molecular profiling and {MRD} detection in Diffuse Large B Cell Lymphoma}, 
   author   = {Krupka, Joanna A and Moutsopoulos, Ilias and Cutmore, Natasha H ... and  Hodson, Daniel J}, 
   journal  = {{medRxiv}}, 
   doi      = {10.1101/2025.04.14.25325806}, 
   pages    = {2025.04.14.25325806}
   }
```

### 🧙‍♂️ About the Name: Aulë

**Aulë the Smith** is one of the Valar (gods) in J.R.R. Tolkien's Middle-earth legendarium. Known as the master craftsman and smith, Aulë created the Dwarves and was renowned for his skill in working with metals, stone, and gems. He was particularly associated with **crafting, precision, and creating tools** that would endure through time.

*"But the delight and pride of Aulë is in the deed of making, and in the thing made, and neither in possession nor in his own mastery"*  
— The Silmarillion

## 💪 Contributors 

[Joanna A. Krupka](https://github.com/ashakru)  
[Ilias Moutsopoulos](https://github.com/emouts)  
[Natasha H Cutmore](https://github.com/nuttynutmore)   
[Daniel J Hodson](https://www.stemcells.cam.ac.uk/people/pi/hodson)  
