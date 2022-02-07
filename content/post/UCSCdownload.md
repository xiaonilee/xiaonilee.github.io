---
title: "Download mRNA lncRNA and miRNA and match clinical information with R from UCSC Xena"
date: 2021-03-12
lastmod: 2021-03-12
draft: false
tags: ["R", "Bioinformatics", "UCSC", "Xena", "TCGA", "Ensembl"]
categories: ["R", "Bioinformatics", "UCSC", "Xena", "TCGA", "Ensembl"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

It's very convenient to get datasets from UCSC Xena. Normally, datasets will need to be processed and transformed with R language. 

Herein, I will recover the protocol of datasets processed.

<!--more-->

### Methods

- **Dataset sources**: [GDC TCGA Head and Neck Cancer(HNSC)](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Head%20and%20Neck%20Cancer%20(HNSC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

- **Object address**: [gene expression RNAseq](https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-HNSC.htseq_counts.tsv.gz), [phenotype](https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-HNSC.GDC_phenotype.tsv.gz) and [survival data](https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-HNSC.survival.tsv.gz)

  ![fig1](fig1.png)

### Prerequisites

```r
library(openxlsx)
library(tidyverse)
library(limma)
library(readr)
```

### Annotation of gene

- **Load data and annotation information**

```r
TCGA_rawdata <- read_tsv("/Users/xiaonili/Downloads/TCGA-HNSC.htseq_counts.tsv.gz")

dim(TCGA_rawdata)

probeMap <- read.table("/Users/xiaonili/Downloads/gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T)
probeMap[1:4,1:4]
```

- **Output**

  ![fig2](fig2.png)

- **ID Reverse**
  
```r
TCGA_gset <- TCGA_rawdata %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  select(gene, starts_with("TCGA") )

TCGA_gset[1:4,1:4]
```

- **Output**
  
  ![fig3](fig3.png)

- **Average replicate genes**

```r
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene))
colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)
write.csv(TCGA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_Countdata_log2+1.csv")
TCGA_gset[1:4,1:4]
```

- **Output**

  ![fig4](fig4.png)

- **Group by patient.id**

```r
TCGA_group_list <- ifelse(as.numeric(substring(colnames(TCGA_gset),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)
```

- **Output**

  ![fig5](fig5.png)

### Recognize mRNA lncRNA and miRNA

- **Load data**

```r
mRNA_info <- read.xlsx("/Users/xiaonili/Downloads/Gene_info.xlsx",sheet = "mRNA_info")
lncRNA_info <- read.xlsx("/Users/xiaonili/Downloads/Gene_info.xlsx",sheet = "lncRNA_info")
miRNA_info <- read.xlsx("/Users/xiaonili/Downloads/Gene_info.xlsx",sheet = "miRNA_info")
```

- **Get geneset for mRNA miRNA and lncRNA** 

```r
## Get data.matrix for mRNA
mRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% mRNA_info$gene_name,]
dim(mRNA_gset)

write.csv(mRNA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_mRNA.csv",quote = F,row.names = T)

## Get data.matrix for lncRNA
lncRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% lncRNA_info$gene_name,]
dim(lncRNA_gset)

write.csv(lncRNA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_lncRNA.csv",quote = F,row.names = T)

## Get data.matrix for miRNA
miRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% miRNA_info$gene_name,]
dim(miRNA_gset)

write.csv(miRNA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_miRNA.csv",quote = F,row.names = T)
```

### Match clinical and survival information with expression

- **Load clinical data**

```r
Phenodata <- read_tsv("/Users/xiaonili/Downloads/TCGA-HNSC.GDC_phenotype.tsv.gz")

Phenodata[1:4,1:4]

Phenodata$submitter_id.samples <- substring(Phenodata$submitter_id.samples,1,15) %>% 
  gsub("-",".",.)
Phenodata[1:4,1:4]
```

- **Output**

  ![fig6](fig6.png)

- **Load survival data**

```r
Sur_data <- read_tsv("/Users/xiaonili/Downloads/TCGA-HNSC.survival.tsv.gz")

Sur_data$sample <- substring(Sur_data$sample,1,15) %>% gsub("-",".",.)
Sur_data[1:4,1:4]
```

- **Output**
  
  ![fig7](fig7.png)

- **Merge data and choose interested col**

```r
Phen_surv <- Phenodata %>%
  inner_join(Sur_data,by = c("submitter_id.samples" = "sample")) %>%
  select(submitter_id.samples,age_at_index.demographic,gender.demographic,
         tumor_grade.diagnoses,neoplasm_histologic_grade,tumor_stage.diagnoses,OS,OS.time)
head(Phen_surv)
```

- **Output**

  ![fig8](fig8.png)

- **match expression with phenodata and do order**

```r
Phen_surv = Phen_surv[match(colnames(TCGA_gset),Phen_surv$submitter_id.samples),]
identical(Phen_surv$submitter_id.samples,colnames(TCGA_gset))
```

- **group**

```r
Phen_surv$group <- TCGA_group_list
Phen_surv = dplyr::select(Phen_surv,submitter_id.samples,group,everything())
write.csv(Phen_surv,"/Users/xiaonili/Downloads/TCGA_HNSC_phenotype.csv")
head(Phen_surv)
```

- **Output**

  ![fig9](fig9.png)

Attach is the [script](UCSCdownload.R).
