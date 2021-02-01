---
title: "Comprehensive Bioinformatics Analysis of Tumor-Infiltrating Immune Cells"
date: 2021-01-31
lastmod: 2021-02-01
draft: false
tags: ["Tumor Immune", "Cancer Research", "Database", "TIMER", "TCIA", "ImmuneCellAI", "EPIC", "xCell"]
categories: ["Cancer Research", "Database"]
author: "Xiaoni"

weight: 1


mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Recent clinical successes of cancer immunotherapy necessitate the investigation of the interaction between malignant cells and the host immune system. 

Herein, I summarized five bioinformatics tools for comprehensive analysis of tumor-infiltrating immune cells.

<!--more-->

### 1. TIMER

- [TIMER](https://cistrome.shinyapps.io/timer/), A web server for comprehensive analysis of tumor-infiltrating immune cells[^1].
- Levels of **six** tumor-infiltrating immune subsets are **pre-calculated** for **10,897 tumors** from **32 cancer types**. 
- 6 major analytic modules that allow users to ***interactively*** explore the associations between immune infiltrates and a wide-spectrum of factors, including `gene expression`, `clinical outcomes`, `somatic mutations`, and `somatic copy number alterations`. 

  ![timer1](timer1.png)

### 2. The Cancer Immunome Atlas

- [TCIA](https://tcia.at/home), provides results of comprehensive immunogenomic analyses of next generation sequencing data (NGS) data for 20 solid cancers from The Cancer Genome Atlas (TCGA) and other datasources[^2].
- Cancer genotypes determine tumor immunophenotypes and tumor escape mechanisms
  ![tcia](tcia2.png)

### 3. ImmuneCellAI

- [ImmuneCellAI](http://bioinfo.life.hust.edu.cn/ImmuCellAI#!/), is a tool to estimate the abundance of 24 immune cells from gene expression dataset including RNA-Seq and microarray data, in which the 24 immune cells are comprised of 18 T-cell subtypes and 6 other immune cells: B cell, NK cell, Monocyte cell, Macrophage cell, Neutrophil cell and DC cell[^3].

  ![immunecellai](immunecellai3.png)

### 4. EPIC

- [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/), is designed to Estimate the Proportion of Immune and Cancer cells from bulk tumor gene expression data[^4]. 
  
  ![epic](epic4.png)

### 5. xCell

- [xCell](https://xcell.ucsf.edu/), is a webtool that performs cell type enrichment analysis from gene expression data for 64 immune and stroma cell types. xCell is a gene signatures-based method learned from thousands of pure cell types from various sources[^5].
  
  ![xCell](xcell5.png)

- [xCellView](http://comphealth.ucsf.edu/xCellView/), visualize cell types enrichments.

  ![xcellview](xcellview52.png)

- For single-cell RNA-seq, to github and use [SingleR](https://github.com/dviraran/SingleR).


[^1]:[Cancer Res. 2017 Nov 1; 77(21): e108–e110.doi: 10.1158/0008-5472.CAN-17-0307](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6042652/)
[^2]:[Cell Rep. 2017 Jan 3;18(1):248-262. doi: 10.1016/j.celrep.2016.12.019](https://pubmed.ncbi.nlm.nih.gov/28052254/)
[^3]:[Adv Sci (Weinh). 2020 Feb 11;7(7):1902880. doi: 10.1002/advs.201902880](https://pubmed.ncbi.nlm.nih.gov/32274301/)
[^4]:[Elife. 2017 Nov 13;6:e26476. doi: 10.7554/eLife.26476](https://pubmed.ncbi.nlm.nih.gov/29130882/)
[^5]:[Genome Biol. 2017 Nov 15;18(1):220. doi: 10.1186/s13059-017-1349-1](https://pubmed.ncbi.nlm.nih.gov/29141660/)
