---
title: "TissGDB: tissue-specific gene database in cancer"
date: 2020-09-10
lastmod: 2020-09-15
draft: false
tags: ["TissGDB", "Cancer", "Database"]
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

TissGDB is the tissue-specific gene annotation database in cancer, aiming to provide a resource or reference for cancer and the related disease studies in the context of tissue specificity.

<!--more-->

## In Brief

- The TissGDB `WEBSITE` is [here](https://bioinfo.uth.edu/TissGDB/).

- **Publication Link** [*Nucleic Acids Res.  2018 Jan 4;46(D1):D1031-D1038. doi: 10.1093/nar/gkx850*](https://pubmed.ncbi.nlm.nih.gov/29036590/)

## Representative Tissue-Specific Gene Expression Resources

- The Human Protein Atlas (`HPA`)
- Tissue-specific Gene Expression and Regulation (`TiGER`)
- Genotype-Tissue Expression (`GTEx`)

## Develop Feature

- `2461` Tissue Specific Genes (TissGenes)
- Across `22` Tissue Types
- Match the `28` Cancer types of The Cancer Genome Atlas (TCGA)

## Analyses

- **Identified** hundreds of TissGenes
- universally **kept** or **lost** tissue-specific gene expression
- cancer type-specific **isoform** expression
- **fusion** with oncogenes or tumor suppressor genes
- **markers** for protective or risk prognosis

## Seven Categories of Annotations

1. TissGeneSummary
2. TissGeneExp
3. TissGene-miRNA
4. TissGeneMut
    - TissGeneSNV
    - TissGeneCNV
    - TissGeneFusions
5. TissGeneNet
     - CePIN
        - Cytoscape
        - R
6. TissGeneProg
     - overall survival, OS
     - relapse free survival, RFS
7. TissGeneClin
     - TissGeneDrug
     - TissGeneDisease

[![TissGDB logo](TissGDB.png)](https://bioinfo.uth.edu/TissGDB/)
