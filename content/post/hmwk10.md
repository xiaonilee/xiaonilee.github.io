---
title: "R Notebook : Application Exercises"
date: 2020-10-15
lastmod: 2020-10-15
draft: false
tags: ["R", "Bioinformatics", Plots]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Bioinformatics Exercises of Application with R

<!--more-->

## In summary

1.strsplit()

2.lapply() row 1, column 2

3.library(stringr) split character

4.merge()

5.match()

6.library(dplyr) delete column

7.library(org.Hs.eg.db) annotation id reverse

8.library(hgu133a.db) annotation id reverse

9.?pData group info

10.library(ggstatsplot)

11.?survfit

12.library(survival)

13.library(survminer)

14.?survfit

15.?ggsurvplot

16.GEOquery download dataset

17.pheatmap::pheatmap

18.?cor

19.'%in%'

20.library(limma) DEG

21.exprs()

### Exercise1. id reverse by org.Hs.eg.db

```r
a=read.table('ensembl_id2symbol.txt')

library(org.Hs.eg.db)
g2s=toTable(org.Hs.egSYMBOL)
g2e=toTable(org.Hs.egENSEMBL)

library(stringr)
a$ensembl_id=str_split(a$V1,'[.]',simplify=T)[,1]

b=merge(a,g2e,by='ensembl_id',all.x=T)
d=merge(b,g2s,by='gene_id',all.x=T)

d=d[match(a$V1,d$V1),]

library(dplyr)
res=select(d,-V1)
```

### Exercise2. id reverse by hgu133a.db

```r
a2=read.table('probe_id2symbol.txt')

# BiocManager::install("hgu133a.db")
library(hgu133a.db)
ls("package:hgu133a.db")

id2s=toTable(hgu133aSYMBOL)
colnames(id2s)
colnames(a2)='probe_id'

b2=merge(a2,id2s,by='probe_id')
res2=b2
```

### Exercise3. get gene expression by CLL and do boxplot

```r
library(CLL)
data("sCLLex")
exprSet=exprs(sCLLex)
dim(exprSet)
exprSet[1:6,1:6]

pd=pData(sCLLex)
colnames(pd)
table(pd$Disease)

library(hgu95av2.db)
ls('package:hgu95av2.db')
id2s=toTable(hgu95av2SYMBOL)
filter(id2s,symbol=='FSCN1')

exprSet['39070_at',]

pd[pd$Disease=='progres.',]

library(ggplot2)
boxplot(exprSet['39070_at',] ~ pd$Disease,col='brown',
        main='fscn1 in CLL',
        xlab='Disease', ylab='expression of fscn1')
ggsave('fscn1boxplot.png')
```

![fscn1boxplot](fscn1boxplot.png)

### Exercise4. get data via TCGA and r analysis

#### setup dataset from TCGA and download

Breast--> TCGA PanCancer Atlas-->FSCN1-->Plots

```r
a=read.table('fscn1.txt',sep = '\t', header = T)
colnames(a)
colnames(a)=c('id','subtype','expression','mut','CNA')

library(ggstatsplot)
ggbetweenstats(data=a, x = subtype, y = expression)
ggsave('fscn1plot.png')
```

![fscn1plot](fscn1plot.png)

### Exercise5. gene in clinical survival

#### get data from oncoLnc.org

- oncoLnc-->TP53-->BRCA-->LOWER/UPPER Percentile of 50-->Submit

```r
a=read.table('BRCA_7157_50_50.csv',sep = ',', header = T)
head(a)
colnames(a)
table(a$Status)

library(ggplot2)
library(survival)
library(survminer)

a$Status=ifelse(a$Status=='Dead',1,0)
fit=survfit(Surv(Days,Status) ~ Group, data = a)
str(fit)
summary(fit)

ggsurvplot(fit, conf.int = F, pval = TRUE,
           linetype = 5,
           ggtheme = theme_survminer())
ggsave('tp53_expression_Survival.png')
```

![tp53_expression_Survival](tp53_expression_Survival.png)

- considering the poor R value, re-analysis with subtype

```r
b=read.table('fscn1.txt',sep = '\t', header = T)
colnames(b)=c('id','subtype','expression','mut','CNA')
colnames(a)

library(stringr)
b$Patient =str_split(b$id,'-01',simplify = T)[,1]

d=merge(a, b, by='Patient', all.x=T)
fit2=survfit(Surv(Days,Status) ~ subtype, data = d)
str(fit2)
summary(fit2)

ggsurvplot(fit2, conf.int = F, pval = TRUE,
           linetype = 5,
           ggtheme = theme_survminer())
ggsave('tp53_Subtype_Survival.png')
```

![tp53_Subtype_Survival](tp53_Subtype_Survival.png)

### Exercise6. Download GSE17215 and plot heatmap
