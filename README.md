# Howell EDA Example

**Author:** Mike Martinez M.S.

**Date Started:** January 27th, 2026

Exploratory data analysis and RGenEDA tutorial for Compute and Conquer February 2026 blog post.



![Status](https://img.shields.io/badge/status-complete-green?style=for-the-badge)



---

## Table of Contents  
- [Background](#background)
- [Software](#software)
- [Data](#data)
- [Methods](#methods)
  - [Data tidying](#data-tidying)
  - [Exploratory data analysis tutorial](#exploratory-data-analysis-tutorial)

## Background

This repository houses a mock exploratory data analysis (EDA) to showcase the usage of [`RGenEDA`](https://github.com/mikemartinez99/RGenEDA) from the publication: [DNA Methylation and Transcription Patterns in Intestinal Epithelial Cells from Pediatric Patients with Inflammatory Bowel Diseases Differentiate Subtypes and Associate With Outcome](https://europepmc.org/article/MED/29031501) (Howell et. al., 2017)
This data consists of 236 intestinal epithelial biopsies from the terminal ileum, ascending, and sigmoid colons from childen with diagnosed IBD (Crohn's (CD): n = 43, ulcerative colitis (UC): n = 23, and no IBD (Control): n = 30). 

Authors of this publication show in their gene expression data a gut segment-specific effect between IBD vs control patients. Based on this *a priori* knowledge, exploratory data analaysis alone should be able to recapitulate this effect before any downstream differential expression modeling occurs. 

This is **not** meant to be a comprehensive re-analysis, rather a tutorial to accompany the February 2026 *Compute and Conquere* blog post. 

---

## Software

All software and R packages are managed by Pixi in the [envs folder](https://github.com/mikemartinez99/Howell-EDA-Example/tree/main/envs). The environment can be reproduced using the `pixi.toml` file using `pixi install`. 

For instructions on installing and configuring Pixi for your system, see [here](https://pixi.prefix.dev/latest/installation/)

---

## Data

Preprocessed data in the form of a count matrix and metadata was obtained from the Gene Expression Atlas under accession [ENA: ERP106487, E-MTAB-5464](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5464/Downloads)

## Methods

### Data tidying

| Script | Purpose |
|--------|---------|
| [`01-Clean-metadata.R`](https://github.com/mikemartinez99/Howell-EDA-Example/blob/main/code/01-Clean-metadata.R) | Tidy metadata for downstream use |

### Exploratory data analysis tutorial

| Script | Purpose |
|--------|---------|
| [`02-EDA.R`](https://github.com/mikemartinez99/Howell-EDA-Example/blob/main/code/02-EDA.R) | Run `RGenEDA` on Howell data to demonstrate EDA |


## Citations

**Howell KJ, Kraiczy J, Nayak KM, et al.** DNA Methylation and Transcription Patterns in Intestinal Epithelial Cells From Pediatric Patients With Inflammatory Bowel Diseases Differentiate Disease Subtypes and Associate With Outcome. Gastroenterology. 2018 Feb;154(3):585-598. DOI: 10.1053/j.gastro.2017.10.007. PMID: 29031501; PMCID: PMC6381389.

**Martinez M.** RGenEDA: A framework for genomic exploratory data analysis. protocols.io; 2025. https://doi.org/10.17504/protocols.io.bp2l6z6rdgqe/v1





