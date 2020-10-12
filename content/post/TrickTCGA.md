---
title: "R re-analysis in TCGA"
date: 2020-10-12
lastmod: 2020-10-12
draft: false
tags: ["R", "Bioinformatics", "TCGA"]
categories: ["Cancer Research", "Database"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Data setup via TCGA.

Re-analysis by R.

<!--more-->

## In Brief

- Setup and download interested data via TCGA database.

- Re-analysis by R.

### Data in TCGA database

- Choose on interested [dataset](https://www.cbioportal.org/).
  
  ![brca](brca.png)

- Select interested data.
  
  ![fscn1](fscn1.png)

- Click on `Plots` and setup
  
  ![plots](plots.png)

- Download [data](plot.txt)
  
### Re-analysis by [R](r.Rmd)

- Read data

```r
a= read.table('plot.txt',sep = '\t', header = T)
colnames(a) = c('id', 'subtype', 'expression', 'mutations', 'CNA')
```

- Plot by R
  
```r
dat = a
library(ggstatsplot)
ggbetweenstats(dat, x = subtype, y = expression)
ggsave('subtypeExp.png')
```

![subtypeExp.png](subtypeExp.png)

- Try other analysis
  
```r
library(ggplot2)
ggplot(dat, aes(x=subtype, y=expression)) +
  geom_violin()

# variance
aov(expression ~ subtype, data=dat)

summary(aov(expression ~ subtype, data=dat))
# output
#              Df Sum Sq Mean Sq F value Pr(>F)
# subtype       4  289.5   72.38   115.5 <2e-16 ***
# Residuals   514  322.0    0.63
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


TukeyHSD(aov(expression ~ subtype, data=dat))
# output
# Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = expression ~ subtype, data = dat)
#
# $subtype
#                                  diff        lwr         upr     p adj
# HER2-enriched-Basal-like  -1.58328100 -1.9422245 -1.22433751 0.0000000
# Luminal A-Basal-like      -1.82596364 -2.0873307 -1.56459662 0.0000000
# Luminal B-Basal-like      -2.09467315 -2.3870044 -1.80234192 0.0000000
# Normal-like-Basal-like    -1.61344005 -2.4101234 -0.81675673 0.0000005
# Luminal A-HER2-enriched   -0.24268264 -0.5610358  0.07567051 0.2273188
# Luminal B-HER2-enriched   -0.51139215 -0.8556211 -0.16716323 0.0005268
# Normal-like-HER2-enriched -0.03015905 -0.8473128  0.78699474 0.9999764
# Luminal B-Luminal A       -0.26870951 -0.5094705 -0.02794855 0.0199389
# Normal-like-Luminal A      0.21252359 -0.5667149  0.99176206 0.9453257
# Normal-like-Luminal B      0.48123310 -0.3089298  1.27139601 0.4552733
```
