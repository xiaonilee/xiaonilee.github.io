---
title: "Example7"
date: 2020-09-10T16:01:23+08:00
lastmod: 2020-09-11T16:01:23+08:00
draft: false
tags: ["UCSC", "Xena", "Breast", "Cancer"]
categories: ["Example", "Bioinformatics", "Database"]
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

## Metastatic breast cancer and sepsis: monocytes

**OUTLINE**
Analysis of the total population of monocytes from patients with metastatic breast cancer (MBC), sepsis, or tuberculosis.

- click on 'VIEW MY DATA' and upload data: 'genomicGDS5819_gene-level.tsv' and 'phenotypeGDS5819.tsv'
  ![uploadedData](./uploadedData.png)
- Back to 'VISUALIZATION'
- Search for a study, choose 'Metastatic Breast Cancer and Sepsis'
- Add phenotype, choose 'Disease'
- Add genes: HLA-DRB6 HLA-DRB4 CD86 IL1B CD83 TNF
- Dataset: My Computer Hub
- choose 'genomicGDS5819_gene-level.tsv'
- Add column, Add genes: RNASE3 NKG7 OLFM4, CHI3L1, CAMP, LOC653600 PGLYRP1 CTSG, DEFA4, MPO, CEACAM6, BPI, AZU1 CYP4F3, ANXA3, MS4A3, COL17A1, MMP9, ARG1, S100P IL18RAP, CA4, Rgr, TNFAIP6
  ![Raw settings](./RawSetting.png)
- for column C and column D, to click on 'Display', Adjust Display Settings, color transformation is be replaced with 'center by column mean: x - column average'
  ![update settings](updateSetting.png)
- figure for the relationship between Disease and column **C** Genes
  ![figure 1](./Disease_ColumnC.png)
- **figure** for the relationship between Disease and column **D** Genes
  ![figure 2](./Disease_ColumnD.png)
- Analysis figure for the relationship between Disease and column **C** Genes, by removing parameter of tuberculosis
  ![data interpretation](./dataInterpretation.png)
- Analysis figure for the relationship between Disease and column **D** Genes, by removing parameter of tuberculosis
  ![data analysis](./dataInterpretation2.png)

**IN SUMMARY**, results provide insight into molecular similarities between monocytes from MBC patients and reprogrammed immunosuppressive monocytes from sepsis patients.
