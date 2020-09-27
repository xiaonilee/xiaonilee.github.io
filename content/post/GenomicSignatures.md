---
title: "UCSC Xena Features: Genomic Signatures"
date: 2020-09-14
lastmod: 2020-09-15
draft: false
tags: ["UCSC", "Xena", "Features", "signatures", "Database"]
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

Genomic signatures, sometimes expressed as a weighted sum of genes, are an algebra over genes.

<!--more-->

## In Brief

- Genomic signatures:
  
  - a **list** of genes
  
  - **weighted** gene(s)

  - an **algebra** over genes

- Example of interpretation
  
  - data from a **publication**[^footnote].
  
  - **gene signature**
    - =FN1 + VIM + ZEB1 + ZEB2 + TWIST1 + TWIST2 + SNAI1 + SNAI2 + CDH2 - CLDN4 - CLDN7 - TJP3 - MUC1 - CDH1

## Protocol

- Click on `VIEW MY DATA` and `Select Data File`.
  
- Choose local file `CMS_subtypes.txt`.
  
- Choose `Phenotypic data: categories or non-genomic in a rectangle (e.g. age, mutation status: 'wt' or 'mutant')`.

- Choose `The first column is sample IDs`.

- Choose `There is other public data in Xena on these samples (e.g. TCGA) and want to connect to it.`.

- Choose Search for study `TCGA Colon and Rectal Cancer (COADREAD)`.

- Click on `IMPORT` and `FINISH`.
  
  ![step 1](step1.png)

- Back to `VISUALIZATION`.

- Click `Search Phenotype` to add Phenotype `CMS_subtype`.

- **Add a new column** and choose `gene expression`: `=FN1 + VIM + ZEB1 + ZEB2 + TWIST1 + TWIST2 + SNAI1 + SNAI2 + CDH2 - CLDN4 - CLDN7 - TJP3 - MUC1 - CDH1`.

- Final input.

  ![step2](step2.png)

- Choose `Chart & Statistics` of Column C `or` Column D to generate figures.
  
  ![figure](figure.png)

**In summary**, the present study proposed that taxonomy of four consensus molecular subtypes (CMS) of colorectal cancer(CRC) reflected significant biological differences in the gene expression-based molecular subtypes.

[^footnote]: [Nat Med. 2015 Nov;21(11):1350-6. doi: 10.1038/nm.3967.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4636487/)
