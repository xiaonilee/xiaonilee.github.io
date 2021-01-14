---
title: "Identification of Therapeutic Targets and Prognostic Biomarkers Among CXC Chemokines in the Renal Cell Carcinoma Microenvironment"
date: 2021-01-05
lastmod: 2021-01-14
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
<hr />

  ![fig24](fig2/fig24.png)
<hr />

  ![fig25](fig2/fig25.png)

- Merge all CXC chemokines together with Google Slides
  - The transcriptional levels of CXCL12(G) were significantly reduced. 
  
  ![figure2](fig2/figure2.png)
### Figure 3. The relative level of CXC chemokines in RCC
- Open [GEPIA](http://gepia.cancer-pku.cn/index.html)
  ![fig31](fig3/fig31.png)

- Click on the button `Multiple Gene Analysis` and `Multiple Gene Comparison`

- Set up the parameters and click the button `Plot`
  ![fig32](fig32/fig32.png)

- Save and edit plot [in cloud](https://plot.ly/create/)
  ![fig33](fig3/fig33.png)

- Click on `Heatmap` to generate the final heatmap
  ![fig34](fig3/fig34.png)
<hr />

- The relative expression of **CXCL14** was the highest.   
  
  ![figure3](fig3/figure3.png)

### Figure 4. Correlation between different expressed CXC chemokines and the pathological stage of RCC patients (GEPIA)
- Open [GEPIA](http://gepia.cancer-pku.cn/index.html)

- Click on `Expression DIY`, choose `Stage plot`
  ![fig41](fig4/fig41.png)

- Set up the parameters for CXCL1 gene and click the button `Plot`  

- Generate result of CXCL1 gene.
  ![fig42](fig4/fig42.png)

- Merge all results of CXC chemokines together with Google Slides.
  ![figure4](fig4/figure4.png)
### Figure 5. The prognostic value of different expressed CXC chemokines in RCC patients in the disease free survival curve (GEPIA)
- Open [GEPIA](http://gepia.cancer-pku.cn/index.html)

- Click `Survival` and choose `Survival Plots`
  ![fig51](fig5/fig51.png)

- Set up the parameters for CXCL1 gene and click the button `Plot`
  ![fig52](fig5/fig52.png) 

- Generate result of CXCL1 gene
  ![fig531](fig5/fig531.png)

- Merge all results of CXC chemokines together with Google Slides.
  ![figure5](fig5/figure5.png)

### Figure 6. The prognostic value of CXC chemokines in RCC patients in the overall survival curve (GEPIA)
- Open [GEPIA](http://gepia.cancer-pku.cn/index.html)

- Click `Survival` and choose `Survival Plots`
  ![fig61](fig6/fig61.png)

- Set up the parameters for CXCL1 gene and click the button `Plot`
  ![fig62](fig6/fig62.png) 

- Generate result of CXCL1 gene
  ![fig631](fig6/fig631.png)

- Merge all results of CXC chemokines together with Google Slides.
  ![figure6](fig6/figure6.png)

### Figure 7. Genetic alteration, neighbor gene network, and interaction analyses of different expressed CXC chemokines in RCC patients
- Open [cBioportal](https://www.cbioportal.org/)
  ![fig71](fig7/fig71.png)

- Set up the parameters and click the button `Submit Query`
  ![fig72](fig7/fig72.png)

- Generate result for Figure 7A.
  ![fig73](fig7/fig73.png)

- Download the **mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM)** and [save](fig7/fig7222.txt).
  ![fig732](fig7/fig732.png)

- To get the correlation heat map of different expressed CXC chemokines in RCC with [R.Script](fig7/fig7B.R).
  ![fig7B](fig7/fig7B.png)

- Opne [STRING](https://string-db.org/)

- Set up the parameters and click on the button `SEARCH`.
  ![fig74](fig7/fig74.png)

- Continue, click on the button `CONTINUE>>`
  ![fig75](fig7/fig75.png)

- Click on the button `Exports`
  ![fig76](fig7/fig76.png)

- Click on the button `Exports` and download the result of Fig7C.
  ![fig7C](fig7/fig7C.png)

- Open [GeneMANIA](http://genemania.org/)
  ![fig77](fig7/fig77.png)

- Enter gene list and click search icon
  ![fig78](fig7/fig78.png)

- Click on the style icon and then click on the save icon to get the result of Fig7D.
  ![fig7D](fig7/fig7D.jpg)

### Table 2. Key regulated factor of CXC chemokines in RCC (TRRUST)
- Open [TRRUST](https://www.grnpedia.org/trrust/)
  ![fig1](tab2/fig1.png)

- Entet query genes and setup parameters, and then click on the button `Submit`
  ![fig2](tab2/fig2.png)

- Generate the result
  ![fig3](tab2/fig3.png)
 
- Click on the relative number to check the details
  ![fig4](tab2/fig4.png)

- Edit the result and get the result of table2
  ![table2](tab2/table2.png)
### Figure 8. The enrichment analysis of different expressed CXC chemokines and 50 most frequently altered neighboring genes in RCC (David 6.8)
- Open [DAVID](https://david.ncifcrf.gov/summary.jsp)
  ![fig1](fig8/fig1.png)

- Follow the protocol step by step.
  ![fig2](fig8/fig2.png)
  ![fig3](fig8/fig3.png)
  ![fig4](fig8/fig4.png)
  ![fig5](fig8/fig5.png)
  ![fig6](fig8/fig6.png)

- Generate the result of bar plot with [R.Script](fig8/fig8.R)
  ![figure8](fig8/figure8.png)

- Alternate, generate the result of dot plot.
  ![figure82](fig8/figure82.png)

### Figure 9. The correlation between different expressed CXC chemokines and immune cell infiltration (TIMER)

- Open [TIMER](https://cistrome.shinyapps.io/timer/)
  ![fig1](fig9/fig1.png)

- Click on `gene`

- Set up the parameters for CXCL1 gene and click the button `Submit`
  ![fig2](fig9/fig2.png)

- get the result of CXCL1 gene
  ![fig310](fig9/fig310.png)
<hr>
  
  ![fig31](fig9/fig31.jpg)

- Merge all results of CXC chemokines together with Google Slides
  
  ![fig99](fig9/fig99.png)

### Table 4. The cox proportional hazard model of CXC chemokines and six tumor-infiltrating immune cells in RCC (TIMER)
- Setup the cox proportional hazard of CXC chemokines
  ![fig1](tab4/fig1.png)

### Supplementary Figure 1. The enrichment analysis of different expressed CXC chemokines and 50 most frequently altered neighboring genes in RCC (Metascape)
- Open [Metascape](https://metascape.org/gp/index.html#/main/step1)
  ![fig1](supfig1/fig1.png)

- Click on the button `Express Analysis`

- Click on the button `Analysis Report Page`
  ![fig2](supfig1/fig2.png)

- Generate the results of Supplementary Figure 1
  ![supfigure1](supfig1/supfigure1.png)
