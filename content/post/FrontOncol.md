---
title: "Identification of Therapeutic Targets and Prognostic Biomarkers Among CXC Chemokines in the Renal Cell Carcinoma Microenvironment"
date: 2021-01-05
lastmod: 2021-01-05
draft: false
tags: ["Recover Paper", "Cancer Research", "Database", "Bioinformatics"]
categories: ["Recover Paper", "Cancer Research", "Database", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Recover Bioinformatics paper w/o code.

<!--more-->

## In Brief
- **Publication** 
  - [*Front Oncol.* 2020 Feb 5. doi: 10.3389/fonc.2019.01555](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7012904/)

- **Keywords**
  - Sixteen CXC chemokines (not including CXCL15)
    - CXCL1, CXCL2, CXCL3, CXCL4, CXCL5, CXCL6, CXCL7, CXCL8, CXCL9, CXCL10, CXCL11, CXCL12, CXCL13, CXCL14, CXCL16, CXCL17 
  - Renal cell carcinoma(RCC), Kidney cancer

- **Methods**
  - `ONCOMINE` 
  - `GEPIA` 
  - `UALCAN` 
  - `cBioPortal` 
  - `GeneMANIA` 
  - `DAVID 6.8` 
  - `Metascape` 
  - `TRRUST` 
  - `LinkedOmics` 
  - `TIMER`

### Figure 1. mRNA levels of CXC chemokines in RCC (ONCOMINE)

- Open [Oncomine database](https://www.oncomine.org/resource/main.html)
  ![fig1](fig1/fig1.png)

- Set up the parameters for CXCL1
  ![fig12](fig1/fig12.png)

- Merge all CXC chemokines together with Google Slides
  - Red: the numbers of datasets with statistically significant mRNA over-expression (red)
  - Blue: downregulated expression of CXC chemokines.
  
  ![figure1](fig1/figure1.png)

### Table 1. The mRNA levels of CXC chemokines in different types of RCC tissues and normal renal tissues at transcriptome level (ONCOMINE)


| TLR    | Type                                           | Fold change | P-value  | t-test | References |
|--------|------------------------------------------------|-------------|----------|--------|------------|
| CXCL3  | `Papillary Renal Cell Carcinoma`               | −2.244      | `1.000`  |`-8.293`| (26)       |
| CXCL6  | Clear Cell Renal Cell Carcinoma                | 30.664      | 4.80E-4  | 4.888  | (27)       |
| CXCL7  | Clear Cell Renal Cell Carcinoma                | −9.410      | `0.991`  | −3.343 | (27)       |
| CXCL9  | Clear Cell Renal Cell Carcinoma                | 2.997       | 1.41E-6  | 7.311  | (28)       |
|        | Clear Cell Renal Cell Carcinoma                | 31.985      | 1.81E-7  | 10.220 | (27)       |
|        | Clear Cell Renal Cell Carcinoma                | 4.648       | 1.26E-5  | 6.703  | (30)       |
|        | Clear Cell Renal Cell Carcinoma                | 7.115       | 2.22E-7  | `7.103`| (29)       |
| CXCL10 | Clear Cell Renal Cell Carcinoma                | 12.873      | 3.10E-12 | 11.075 | (27)       |
|        | Clear Cell Renal Cell Carcinoma                | 5.447       | 5.90E-8  | 9.505  | (30)       |
|        | Hereditary Clear Cell Renal Cell Carcinoma     | 11.612      | 9.94E-11 | 9.867  | (29)       |
|        | Non-Hereditary Clear Cell Renal Cell Carcinoma | 5.897       | 4.41E-7  | 6.000  | (29)       |
| CXCL11 | Hereditary Clear Cell Renal Cell Carcinoma     | 2.994       | 9.61E-9  | 7.199  | (29)       |
|        | Clear Cell Renal Cell Carcinoma                | 6.303       | 1.26E-4  | 5.000  | (30)       |
|        | Clear Cell Renal Cell Carcinoma                | 20.691      | 8.45E-4  | 5.712  | (27)       |
| CXCL13 | Clear Cell Renal Cell Carcinoma                | `9.934`     | `0.002`  | `3.633`| (27)       |
|        | Hereditary Clear Cell Renal Cell Carcinoma     | `1.921`     | `7.84E-4`| `3.447`| (29)       |
| CXCL16 | Clear Cell Renal Cell Carcinoma                | 5.797       | 5.82E-4  | 6.812  | (27)       |
|        | Clear Cell Renal Cell Carcinoma                | 2.212       | 7.86E-4  | 3.932  | (28)       |
### Figure 2. The transcription of CXC chemokines in RCC (UALCAN) 
- Open [UALCAN database](http://ualcan.path.uab.edu/http://ualcan.path.uab.edu/)
  ![fig21](fig2/fig21.png)

- Click on the button `TCGA analysis` and explore the result for CXCL1 gene.
  ![fig22](fig2/fig22.png)

- Click on the button `Expression` and generate the expression result of CXCL1 gene.
  ![fig23](fig2/fig23.png)
  ![fig24](fig2/fig24.png)
  ![fig25](fig2/fig25.png)

- Merge all CXC chemokines together with Google Slides
  ![figure2](fig2/figure2.png)
### Figure 3. The relative level of CXC chemokines in RCC

### Figure 4. Correlation between different expressed CXC chemokines and the pathological stage of RCC patients (GEPIA)

### Figure 5. The prognostic value of different expressed CXC chemokines in RCC patients in the disease free survival curve (GEPIA)

### Figure 6. The prognostic value of CXC chemokines in RCC patients in the overall survival curve (GEPIA)

### Figure 7. Genetic alteration, neighbor gene network, and interaction analyses of different expressed CXC chemokines in RCC patients

### Figure 8. The enrichment analysis of different expressed CXC chemokines and 50 most frequently altered neighboring genes in RCC (David 6.8)

### Figure 9. The correlation between different expressed CXC chemokines and immune cell infiltration (TIMER)

