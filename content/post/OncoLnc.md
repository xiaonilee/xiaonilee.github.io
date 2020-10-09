---
title: "OncoLnc: A straightforward Tool"
date: 2020-10-09
lastmod: 2020-10-09
draft: false
tags: ["R", "Bioinformatics", "TCGA", "OncoLnc"]
categories: ["Cancer Research", "Web Tool"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

OncoLnc Web Tool: Interactively exploring survival correlations, and for downloading clinical data coupled to expression data for mRNAs, miRNAs, or lncRNAs from TCGA.

<!--more-->

## In Brief

- Offical [website](http://www.oncolnc.org/) and [Publication Link](https://peerj.com/articles/cs-67/).
  
- Data Wrangling and Re-analysis with R Package.

### OncoLnc Usage

- Enter Interested Gene, **DONSON**, and click on `Submit`.

![DONSON](DONSON.png)

- Choose Interested cancer, **KIRC**, and click on [**Yes Please!**](http://www.oncolnc.org/kaplan/?cancer=KIRC&gene_id=29980&raw=DONSON&species=mRNA) of the cancer. Herein, DONSON gene ranked in #1.

![cancer](cancer.png)

- Setup value for Lower Percentile and Upper Percentile with the purpose of dividing the patients without overlapping slices.

- Click on [Submit](http://www.oncolnc.org/kaplan/?lower=25&upper=25&cancer=KIRC&gene_id=29980&raw=DONSON&species=mRNA).
  
![plot](plot.png)

- Click on [Click Here](KIRC_29980_25_25.csv), to get the excel file of this data.

### Re-analysis with R
