---
title: "OncoLnc: A straightforward Tool"
date: 2020-10-09
lastmod: 2020-10-16
draft: false
tags: ["R", "Bioinformatics", "TCGA", "OncoLnc", "Plots"]
categories: ["Cancer Research", "Web Tool"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

**OncoLnc Web Tool**: Interactively exploring survival correlations, and for downloading clinical data coupled to expression data for mRNAs, miRNAs, or lncRNAs from ***TCGA***.

<!--more-->

## In Brief

- Offical [website](http://www.oncolnc.org/) and [Publication Link](https://peerj.com/articles/cs-67/).
  
- Data Wrangling and Re-analysis with R Package.

- Summary

### OncoLnc Usage

- Enter Interested Gene, **DONSON**, and click on `Submit`.

![DONSON](DONSON.png)

- Choose Interested cancer, **KIRC**, and click on [**Yes Please!**](http://www.oncolnc.org/kaplan/?cancer=KIRC&gene_id=29980&raw=DONSON&species=mRNA) of the cancer. Herein, DONSON gene ranked in #1.

![cancer](cancer.png)

- Setup value for Lower Percentile and Upper Percentile with the purpose of dividing the patients without overlapping slices.

- Click on [Submit](http://www.oncolnc.org/kaplan/?lower=25&upper=25&cancer=KIRC&gene_id=29980&raw=DONSON&species=mRNA).
  
![plot](plot.png)

- Click on [Click Here](KIRC_DONSON.csv), to get the excel file of this data.

### Re-analysis with R

#### Reading Data and Check Basic Information

```r
a=read.table('KIRC_DONSON.csv', header = T, sep = ',', fill = T)

colnames(a)
output
[1] "Patient"    "Days"       "Status"     "Expression" "Group"  

table(a$Status)
output
Alive  Dead
  151   109

fivenum(a$Expression)
output
[1]  15.780  81.605 155.420 266.310 905.630
```

#### Load Packages/Library and Make Basic Plot

```r
library(ggstatsplot)
library(ggplot2)
library(survival)
#install.packages('survminer')
library(survminer)

# 'Dead':1; 'Alive':0
a$Status=ifelse(a$Status=='Dead',1,0)
```

#### Expression vs Group

```r
ggbetweenstats(data=a, x=Group, y=Expression)
ggsave('ExpressionGroup.png')
```

![plot1](plot1.png)

#### Basic Survival Curve

```r
fit = survfit(Surv(Days,Status)~Group,data = a)
fit
summary(fit)

ggsurvplot(fit,conf.int = F, pval = TRUE)
ggsave('survival_R.png')
```

![survival_R](survival_R.png)

#### Survival Curve With More Information In Figure

```r
#  Better Display
ggsurvplot(fit,palette = c("#E7B800","#2E9FDF"),risk.table = TRUE,pval = TRUE,
conf.int = TRUE,xlab="Time in months",
ggtheme = theme_light(),
ncensor.plot=TRUE)
ggsave('survival_R1.png')
```

![survival_R1](survival_R1.png)

Herein, [R.script](OncoLnc.Rmd) is reused to make survival curves.

### Summary

#### Overall

- The code was run with Python 2.7.5, NumPy 1.7.1, and rpy2 2.5.6.

- It can require upwards of 6GB of RAM.

- OncoLnc runs on Django 1.8.2, Python 2.7, matplotlib 1.2.1, NumPy 1.7.1, rpy2 2.5.6, uses the SQLite3 database engine, and utilizes Bootstrap CSS and JavaScript, and Font Awesome icons.

#### [Reproduce information](https://github.com/OmnesRes/onco_lnc)

- Put the expression files to the correct locations.

- Run the desired Cox regressions.

- Go to the cancer of interest and run the cox_regression.py file from the command line.
