---
title: "FireBrowse Notebook"
date: 2020-09-24
lastmod: 2020-09-24
draft: false
tags: ["TCGA", "firehose"]
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

[FireBrowse Homepage](http://firebrowse.org/) and [*BROAD GDAC*](http://gdac.broadinstitute.org/).

Analysis with FIREHOSE plus TCGA UCSC Xena.

<!--more-->

## In Brief

- FIREHOSE Usage.
  
- FIREHOSE plus TCGA UCSC Xena.

## FIREHOSE Usage

### 1. Interested Disease from FIREHOSE

- 1.1 Choose **Breast invasive carcinoma**.

- 1.2 Click on [***Browse** of Analyses*](http://firebrowse.org/?cohort=BRCA) from Diseases list.

  ![browsedisease1](browsedisease1.png)

### 2. Get Reports of Analysis

- Click on [***BRCA***](http://gdac.broadinstitute.org/runs/analyses__latest/reports/cancer/BRCA/) of `TCGA data version 2016_01_28 for BRCA` to get the reports.
  
  ![getreport2](getreport2.png)

### 3. Result Analysis

- In the **Results** section, choose interested index.
  
  ![resultanalysis1](resultanalysis1.png)

## FIREHOSE plus TCGA UCSC Xena

### Example 1

- Click on [`View Report`](http://gdac.broadinstitute.org/runs/analyses__latest/reports/cancer/BRCA/Correlate_CopyNumber_vs_mRNAseq/nozzle.html) of **Correlations between copy number and mRNAseq expression**.
  
  ![resultanalysis2](resultanalysis2.png)

- Click on **Correlation results**.
  
    ![viewreport1](viewreport1.png)

#### 1. Interested gene, ERBB2

- Focus on `ERBB2` gene, which is very [famous](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ERBB2&keywords=ERBB2) gene.

  ![ERBB2](ERBB2.png)

#### 2. Settings for ERBB2

- Set up ERBB2 gene on UCSC Xena according to firehose.

  ![TCGAset](TCGAset.png)

#### 3. Figure 1

- As show in figure1, `r = 0.8617 (p = 0.000)`, which is almost the same as that of firehose Database.

  ![figure1](figure1.png)

### Example 2

- Click on [`View Report`](http://gdac.broadinstitute.org/runs/analyses__latest/reports/cancer/BRCA/Correlate_Methylation_vs_mRNA/nozzle.html) of **Correlation between mRNA expression and DNA methylation**.
  
  ![resultanalysis3](resultanalysis3.png)

- Click on **Negative Correlation between Methylation and Gene Expression**.
  
  ![viewreport2](viewreport2.png)

#### 1. Interested gene, MKRN3

- Focus on `MKRN3` gene, which ranked No.1 in the list.
  
  ![MKRN3](MKRN3.png)

#### 2. Settings for MKRN3

- Set up MKRN3 gene on UCSC Xena according to firehose.

  ![set2](set2.png)

#### 3. Figure 2

- As show in figure2, `r = -0.8392 (p = 1.538e-232)`, which is close to that of firehose Database.

  ![figure2](figure2.png)

## In Summary

- Different database resources, but similarly conclusion.

- Start with preliminary conclusion get from bioinformatics analysis, we can go through to verify by performacing related experiments.
  
- Start with verified data from experiments, we can furthermore get supported by bioinformatics analysis.
  