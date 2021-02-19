---
title: "Local UCSC Xena Data Hub"
date: 2020-09-10
lastmod: 2021-02-19
draft: false
tags: ["UCSC", "Xena", "Breast", "Cancer", "Database"]
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

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

Load the private data to Local Xena Data Hub.

<!--more-->

## In Brief

The Local Xena Browser is the same as the public browser, which shows same `VISUALIZATION` and `TOOLS` no matter where the data is located.

An Application example: Analysis of the total population of monocytes from patients with metastatic breast cancer (MBC), sepsis, or tuberculosis.

## Installing a [Local Xena Hub](https://ucsc-xena.gitbook.io/project/local-xena-hub/getting-started#installing-a-local-xena-hub)

`Double click` on the download to finish the install.

## Dateset Resource

Download related files of GEO accession: [GSE65517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65517).

## Upload Datasets

- Click on `VIEW MY DATA` and upload `local` data: `genomicGDS5819_gene-level.tsv` and `phenotypeGDS5819.tsv`.

- Name study: `Metastatic Breast Cancer and Sepsis`.
  
  ![dataloaded](dataloaded.png)

## Add Interested Columns Step by Step

- Back to `VISUALIZATION`.

- Search for a study, choose `Metastatic Breast Cancer and Sepsis`.
  
  ![study](study.png)
  
- Add phenotype, choose `Disease`.

- Add genes: `HLA-DRB6`, `HLA-DRB4`, `CD86`, `IL1B`, `CD83`, `TNF`.

- Dataset: My Computer Hub.

- Choose `genomicGDS5819_gene-level.tsv`.

- Add column, Add genes: `RNASE3`, `NKG7`, `OLFM4`, `CHI3L1`, `CAMP`, `LOC653600`, `PGLYRP1`, `CTSG`, `DEFA4`, `MPO`, `CEACAM6`, `BPI`, `AZU1`, `CYP4F3`, `ANXA3`, `MS4A3`, `COL17A1`, `MMP9`, `ARG1`, `S100P` `IL18RAP`, `CA4`, `Rgr`, `TNFAIP6`.
  
  ![rawsetting](raw_set.png)

## Adjust Display

- For column `C` and column `D`, to click on `Display`, on the window of **Adjust Display Settings**, color transformation is be replaced with `center by column mean: x - column average`.
  
  ![adjustdisplaysettings](adjustdisplaysettings.png)

- After adjusting, here is the `NEW` display.
  
  ![updatesetting](update_set.png)

## Output

- **Analysis figure** for the relationship between `Disease` and column `C` Genes.
  
  ![figure1](Disease_ColumnC.png)

- **Analysis figure** for the relationship between `Disease` and column `D` Genes.
  
  ![figure2](Disease_ColumnD.png)

- **Analysis figure** for the relationship between `Disease` and column `C` Genes, by removing parameter of tuberculosis.
  
  ![datainterpretation](dataInterpretation.png)

- **Analysis figure** for the relationship between `Disease` and column `D` Genes, by removing parameter of tuberculosis.

  ![dataanalysis](dataInterpretation2.png)

**IN SUMMARY**, results provide insight into molecular `similarities` between monocytes from MBC patients and reprogrammed immunosuppressive monocytes from sepsis patients[^footnote].

[^footnote]: [PLoS One.2015 doi: 10.1371/journal.pone.0127028.](https://pubmed.ncbi.nlm.nih.gov/25992611/)
