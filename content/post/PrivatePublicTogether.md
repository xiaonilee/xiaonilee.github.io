---
title: "UCSC Xena: Private and Public Data Put Together"
date: 2020-09-11
lastmod: 2020-10-07
draft: false
tags: ["UCSC", "Xena", "Prostate", "Cancer", "Database"]
categories: ["Cancer Research", "Database"]
author: "Xiaoni"

weight: 1

# You can also close(false) or open(true) something for this content.
# P.S. comment can only be closed
# comment: false
# toc: false

# You can also define another contentCopyright. e.g. contentCopyright: "This is another copyright."
contentCopyright: ''
# reward: false
mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Put private data and public Xena Data together via UCSC Xena Hub .

<!--more-->

## In Brief

- Load local data.

- Set up public UCSC Xena data.
  
- Local data is from publication [*Cell Rep. PMID# 29617662*](https://www.sciencedirect.com/science/article/pii/S2211124718303954).

- TCGA database is `TCGA Prostate Cancer (PRAD)`.

## Load TCGA Database

- Study: TCGA Prostate Cancer (PRAD).

- Add 3 columns for `ERG`: `Gene expression`, `Copy number`, `Exon expression`.

- Add column for `TMPRSS2`, `Copy number`.

- Add column for `chr21`, `Copy number`.
  
  ![step 0](step0.png)

## Adding In Data From the Literature: ERG-TMPRSS2 fusion
  
- Click on `View My Data` in the top navigation bar and load data.
  
  ![step 1](step1.png)

- Load as `Phenotype` data.
  
  ![step 2](step2.png)

- `Sample ID` looks like: TCGA-VP-A87D-01.
  
  ![step 3](step3.png)

- Should suggest Study as `TCGA Prostate Cancer (PRAD)`.
  
  ![step 4](step4.png)

- `IMPORT` and `FINISH`.
  
  ![step 5](step5.png)

- Back to `Visualization`.
  
  ![step 6](step6.png)

- Add **Phenotype** `ERG-TMPRSS2 fusion` and generate the final input.
  
  ![step 7](step7.png)

**In summary**, as shown in the last image:

- About half of the samples have high expression of ERG ( Column `B`).

- The high expression of ERG is specifically located in the exons from 3' end of genes ( Column `C`).

- Samples with over-expression of ERG shows a deletion in chromosome 21 ( Column `D`).
  
- Fusion calls agree with expression and Copy number data ( Column `EFG`).
