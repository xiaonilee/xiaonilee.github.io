---
title: "UCSC Xena Hub: Load Data From Public Database"
date: 2020-09-08
lastmod: 2020-10-07
draft: false
tags: ["UCSC", "Xena", "GEO", "Database"]
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

Load and analysis data from publication resources, such as [GEO](https://www.ncbi.nlm.nih.gov/geo/).

<!--more-->

## In Brief

- Download data of GEO.

- Load data step by step.

## Protocol

### Download Data

- Enter and search GEO accession: `GSE121720`.
  
  ![searchgse](searchGSE121720.png)

- Download and save related data for GSE121720.
  
  ![GSE121720](GSE121720.png)

### Load data

- Click on `VIEW MY DATA`.
  
- Click on `LOAD DATA` on the right of **My computer hub**.

  ![step1](step1.png)

- Click on `Select Data File`.

- Choose `GSE121720_phenotype.txt`.

- Click on 'NEXT'

- Choose `Phenotypic data: categories or non-genomic in a rectangle (e.g. age, mutation status: 'wt' or 'mutant')`.

  ![step2](step2.png)

- Click on 'Next'

- Choose `The first row is sample IDs`

  ![step3](step3.png)

- Click on 'Next'

- Choose `These are the first data on these samples`.

- Click on `IMPORT`

  ![step4](step4.png)

- Choose `LOAD MORE DATA`.

  ![step5](step5.png)

- Click on `Select Data File`.

- Choose `GSE121720_RNAseq_expression_matrix_TPM.txt`.

- Click on 'Next'

- Choose `Genomic data: numbers in a rectangle (e.g. expression)`.

  ![step6](step6.png)

- Click on 'Next'

- Choose `The first row is sample IDs`.

- Click on 'Next'

- Choose `I have loaded other data on these samples and want to connect to it`.

- Select a study and choose `My Study`.

- Click on `ADVANCED`

  ![step7](step7.png)

- Click on `IMPORT`

  ![step8](step8.png)

- Click on `FINISH`

  ![step9](step_over.png)

- Click on `VISUALIZE` to continue to do interested analysis.
  
### Output

  ![figure](IDH1Mutation_Expression.png)
  