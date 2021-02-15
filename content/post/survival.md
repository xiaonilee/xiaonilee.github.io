---
title: "Survival Analysis with data from GSE database"
date: 2021-02-15
lastmod: 2021-02-15
draft: false
tags: ["R", "Bioinformatics", "GSE", "Survival analysis"]
categories: ["R", "Bioinformatics", "GSE", "Survival analysis"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Class topics - Survival Analysis.

<!--more-->

## In summary

**Data set**: Any from GSE database, herein, it is `GSE8894`.

**Interested Genes**: Any one(ones) from the data set, herein, they are `Cx43` and `PCDH7`.

**Analysis**: Recurrence free survival

**Result show**: Extended Data Fig.2e from [Nature. 2016 May 26; 533(7604): 493–498](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5021195/)

### Download expression matrix and clinical information

```r
library(GEOquery)
gse8894 = getGEO('GSE8894', destdir = '.', AnnotGPL = F, getGPL = F)
#Check class
class(gse8894)
```


```r
eSet <- gse8894[[1]]
dat_expr <- exprs(eSet)
dim(dat_expr)
head(dat_expr[,1:4])
```


```r
#log2()
dat_expr <- log2(dat_expr+1)
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
colnames(phenoDat)
table(phenoDat$`cell type:ch1`)

#filter only Adenocarcinoma
pd <- phenoDat[phenoDat$`cell type:ch1` =='Adenocarcinoma',]
colnames(pd)
pd = pd[,c(2,40,41)]

# rename colnames
colnames(pd) = c('id','time','event')
meta = pd # survival information
head(meta)

dat_expr = as.data.frame(dat_expr)

exp = dat_expr[,colnames(dat_expr) %in% pd$id] # with "geo_accession" and "probe_id"

GPL=eSet@annotation

# convert gene symbol and probe_id
if(F){
  if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
  library(hgu133plus2.db)
  ls("package:hgu133plus2.db")
  ids <- toTable(hgu133plus2SYMBOL)
  head(ids)
}else if(F){
  getGEO(GPL)
  ids = data.table::fread(paste0(GPL,".soft"),header = T,skip = "ID",data.table = F)
  ids = ids[c("ID","Gene Symbol"),]
  colnames(ids) = c("probe_id")
}else if(T){
  ids = idmap(GPL,type = "bioc")
}

head(ids)
probes_anno=ids
head(probes_anno) 

# Filter expression matrix based on annotation
# and remove the duplicated gene symbols
genes_exp <- filterEM(exp,probes_anno)


## ======================Check up interested genes======================
head(genes_exp)
genes_exp['PCDH7',]
genes_exp['GJA1',]

## Extract and save
mat=genes_exp[c('GJA1','PCDH7'),] #expression information
save(mat,meta,file = 'input_for_survial.Rdata')
list.files()


## ======================Analysis results============================
rm(list = ls())
options(stringsAsFactors = F)
library(survival)
library(survminer)
load(file = 'input_for_survial.Rdata')
identical(colnames(mat),rownames(meta))

e=cbind(meta,t(mat))
head(e)
e$time = as.numeric(e$time)
e$event = as.numeric(e$event)
str(e)
colnames(e)

# Determine the optimal cutpoint for each variable using 'maxstat.'
tmp.cut <- surv_cutpoint(
  e,
  time = "time",
  event = "event",
  variables = c('GJA1','PCDH7')
)
summary(tmp.cut)
plot(tmp.cut, "GJA1", palette = "npg")
# Divide each variable values based on the cutpoint returned by surv_cutpoint().
tmp.cat <- surv_categorize(tmp.cut, labels = c("low","high"))
head(tmp.cat)
table(tmp.cat[,3:4])

kp=tmp.cat$GJA1 == tmp.cat$PCDH7
table(kp)
tmp.cat=tmp.cat[kp,]

#plot
library(survival)
fit <- survfit(Surv(time, event) ~ GJA1+PCDH7,data = tmp.cat)
ggsurvplot(
  fit,
  palette = c("brown", "green"),
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,120),
  break.time.by = 20,
  legend=c(0.5, 0.9),
  legend.title = '',
  font.legend = c(16, "bold.italic", "black"),
  risk.table.y.text.col = T,
  risk.table.y.text = FALSE,
  xlab = "Months",
  ylab = "Metastasis free survival",
  font.x = c(16, "bold", "black"),
  font.y = c(16, "bold", "black"),
  font.title = c(16, "bold", "black"),
  font.subtitle = c(16, "bold", "black"),
  title = "lung adenocarcinoma",
  subtitle = "(GSE8894)"
) 
```

### Output

![fig1](fig1p2.png)
