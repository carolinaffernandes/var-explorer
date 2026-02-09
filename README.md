[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18563349.svg)](https://doi.org/10.5281/zenodo.18563349)


## Overview

This pipeline is an **exploratory framework for the analysis of variant call files (.vcf)** and vcf files annotated with **ANNOVAR**.

It summarizes functional impact, genomic distribution, substitution patterns, and recurrent genes across samples or experimental groups, generating:

- publication-ready plots  
- tabular summaries  
- an integrated HTML report  

---

## Setup

### 1. Configure the analysis

Fill the `config.yaml` file according to the template provided in the repository.

### 2. Create the conda environment

Run:

    conda env create -f environment.yaml
    conda activate var-explorer

### 3. Project structure

- `results/` → analysis outputs  
- `scripts/` → pipeline source code  
- `run_pipeline.sh` → main execution script  

- `scripts/summary.py` → AF and DP summary by sample
- `scripts/report.py` → Variant attribute summarization
---

## Running the pipeline

Execute from the repository root:

```bash
bash run_pipeline.sh
```

You can also run a dry test without executing commands:

```bash
bash run_pipeline.sh --dry-run
```

---

## Variant impact classification

Variants are categorized into four biological impact levels:

| Class | Biological meaning |
|-------|--------------------|
| **HIGH** | Likely protein-disrupting (stopgain, frameshift, splicing) |
| **MODERATE** | Amino acid change (missense) with potential functional effect |
| **LOW** | Synonymous change with minimal functional impact |
| **MODIFIER** | Non-coding or uncertain functional consequence |

This classification is conceptually consistent with **VEP / SnpEff impact annotations** and is commonly used in **somatic variant prioritization**.

---

## Generated plots

The pipeline produces:

- **Genomic region distribution**
- **Exonic functional classification**
- **Variant impact levels**
- **Top recurrent genes**
- **Chromosomal distribution**
- **Base substitution types (SBS-6, COSMIC-style ordering)**

All figures are exported as **publication-ready PNG files** and summarized in an **interactive HTML report**.

---

## Documentation and references

#### Variant annotation tools

If you use this pipeline in scientific work, please also consider citing the annotation frameworks conceptually related to the impact classification:

**SnpEff**  
Cingolani P, Platts A, Wang LL, et al.  
*A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff.*  
Fly (Austin). 2012;6(2):80–92.  
https://doi.org/10.4161/fly.19695

**Ensembl Variant Effect Predictor (VEP)**  
McLaren W, Gil L, Hunt SE, et al.  
*The Ensembl Variant Effect Predictor.*  
Genome Biology. 2016;17:122.  
https://doi.org/10.1186/s13059-016-0974-4

---

### Python libraries

- BeautifulSoup documentation:  
  https://www.crummy.com/software/BeautifulSoup/bs4/doc.ptbr/

- BeautifulSoup official docs:  
  https://beautiful-soup-4.readthedocs.io/en/latest/

- Tutorial:  
  https://dev.to/leapcell/scrape-like-a-pro-beautifulsoup-python-full-tutorial-nj4

