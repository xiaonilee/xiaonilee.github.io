---
title: "R Notebook : Application Exercises"
date: 2020-10-15
lastmod: 2020-10-16
draft: false
tags: ["R", "Bioinformatics", "Plots", "TCGA", "OncoLnc"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Bioinformatics Exercises of Application with [R](hmwk10.Rmd).

<!--more-->

### In summary

1.strsplit()

2.`lapply()` row 1, column 2

3.`library(stringr)` split character

4.merge()

5.match()

6.library(dplyr) delete column

7.library(org.Hs.eg.db) annotation `id reverse`

8.library(hgu133a.db) annotation id reverse

9.?pData `group info`

10.`library(ggstatsplot)`

11.?survfit

12.`library(survival)`

13.library(survminer)

14.?survfit

15.?ggsurvplot

16.`GEOquery` download dataset

17.`pheatmap::pheatmap`

18.`?cor`

19.`'%in%'`

20.library(limma) `DEG`

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

```r
# download GEO dataset
# BiocManager::install('GEOquery')
library(GEOquery)
gse17215 = getGEO('GSE17215', destdir = '.', AnnotGPL = F, getGPL = F)
class(gse17215)
length(gse17215)

a=gse17215[[1]]
dat=exprs(a)

dim(dat)
colnames(dat)
rownames(dat)

library(hgu133a.db)
library()
id2s=toTable(hgu133aSYMBOL)

dat = dat[id2s$probe_id,]
id2s$median=apply(dat,1,median)
id2s=id2s[order(id2s$symbol,id2s$median,decreasing = T),]
id2s=id2s[!duplicated(id2s$symbol),]
dat=dat[id2s$probe_id,]
rownames(dat)=id2s$symbol

# gene list
gl='ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T'

gl = strsplit(gl, ' ')[[1]]

table(gl %in% rownames(dat))

gl=gl[gl %in% rownames(dat)]
dat = dat[gl,]
dat = log2(dat)

pheatmap::pheatmap(dat, scale = 'row')
ggsave('pheatmap_row.png')
```

![pheatmap_row](pheatmap_row.png)

### Exercise7. Download GSE24673 to do correlation and plot heatmap

```r
library(GEOquery)
gse24673=getGEO('GSE24673', destdir = '.', AnnotGPL = F, getGPL = F)

class(gse24673)
length(gse24673)

a = gse24673[[1]]
dat=exprs(a)
dim(dat)
colnames(dat)

pd=pData(a)

# create group_list, according to pd 'source name'.
group_list=c('rbc','rbc','rbc',
             'rbn','rbn','rbn',
             'rbc','rbc','rbc',
             'normal','normal')


M=cor(dat)
colnames(M)
pheatmap::pheatmap(M)
ggsave('pheatmap2.png')

tmp=data.frame(group=group_list)
rownames(tmp)=colnames(M)
pheatmap::pheatmap(M, annotation_col = tmp)
ggsave('pheatmap3.png')
```

![pheatmap3](pheatmap3.png)

### Exercise8. Annotation package install

```r
# ref:https://gist.github.com/xiaonilee/f582c75173122d3625afc34f638b7b30
BiocManager::install('hugene10sttranscriptcluster.db')
```

### Exercise9. Download GSE42872 to compute ave/sd/mad of probes/genes

```r
library(GEOquery)
gse42872 = getGEO('GSE42872', destdir = '.', AnnotGPL = F, getGPL = F)
a = gse42872[[1]]
dat=exprs(a)
pd=pData(a)
colnames(dat)
rownames(dat)

sort(apply(dat, 1, mean), decreasing = T)[1]
# output
# 7978905
# 14.53288

sort(apply(dat,1,sd), decreasing = T)[1]
# output
# 8133876
# 3.166429

sort(apply(dat, 1, mad), decreasing = T)[1]
# output
# 8133876
# 4.268561


# check annotation package https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872
library(hugene10sttranscriptcluster.db)
ls('package:hugene10sttranscriptcluster.db')
id2s=toTable(hugene10sttranscriptclusterSYMBOL)


max_mean = id2s[id2s$probe_id %in% 7978905,]
table(id2s$probe_id %in% 7978905)

max_sd = id2s[id2s$probe_id %in% 8133876,]
max_mad = id2s[id2s$probe_id %in% 8133876,]
```

### Exercise10. Download GSE42872 to do DEG

```r
library(GEOquery)
gse42872 = getGEO('GSE42872', destdir = '.', AnnotGPL = F, getGPL = F)
a = gse42872[[1]]
dat=exprs(a)
dim(a)
pd=pData(a)
colnames(dat)

colnames(pd)
# strsplit(pd$title,' ')
# str_split(pd$title,' ')
# strsplit(pd$title,' ')[[1]]
# strsplit(pd$title,' ')[[1]][4]

group_list=unlist(lapply(pd$title, function(x){
  str_split(x, ' ')[[1]][4]
}))
group_list

# design matrix
library(limma) # Linear Models for Microarray Data
design = model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(dat)

contrast.matric = makeContrasts(paste0(unique(group_list), collapse = '-'), levels = design)

fit = lmFit(dat,design) # Linear Model for Series of Arrays
fit = contrasts.fit(fit,contrast.matric)
fit = eBayes(fit)

tempOutput = topTable(fit,coef = 1, n=Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
```
