---
title: "Chapter 18. Advanced methods for missing data"
date: 2020-11-16
lastmod: 2020-11-16
draft: false
tags: ["R", "R in Action", "Bioinformatics", "Book"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Notebook of Reading Books: R in Action_Chapter 18.

<!--more-->

## This chapter covers

- Identifying missing data

- Visualizing missing data patterns

- Complete-case analysis

- Multiple imputation of missing data

### 18.1. Steps in dealing with missing data

- Identidy missing data.
- Examine the causes of missing data.
- Delete missing data.

  ![fig181](fig181.png)

### 18.2. Identifying missing values

### 18.3. Exploring missing-values patterns

#### 18.3.1. Tabulating missing values

- md.pattern()

#### 18.3.2. Exploring missing data visually

- Figure 18.2. `md.pattern()`-produced plot of missing-values for the sleep dataset in the `VIM` package.
  
  ![fig182](fig182.png)

- Figure 18.2-2. `aggr()`-produced plot of missing-values patterns for the sleep dataset.
  
  ![fig1822](fig182-2.png)

- Figure 18.3. Matrix plot of actual and missing values by case (row) for the sleep dataset. 
  - The matrix is sorted by BodyWgt.

  ![fig1833](fig183.png)

- Figure 18.4. Scatter plot between amount of dream sleep and length of gestation, with information about missing data in the **margins**.
  - with marginplot()

  ![fig184](fig184.png)

#### 18.3.3. Using correlations to explore missing values

### 18.4. Understanding the sources and impact of missing data

### 18.5. Rational approaches for dealing with missing data

### 18.6. Complete-case analysis(listwise deletion)

### 18.7. Multiple imputation

- Figure 18.5. Steps in applying multiple imputation to missing data via the mice approach.

  ![fig185](fig185.png)

### 18.8. Other approaches to missing data

#### 18.8.1. Pairwise deletion

#### 18.8.2. Simple(nonstochastic) imputation

Attach is the [Script](chapter18.R) of chapter18.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 18.2. Identifying missing values
# install.packages("VIM")
data(sleep, package="VIM")

# list the rows that do not have missing values
sleep[complete.cases(sleep),]

# list the rows that have one or more missing values
sleep[!complete.cases(sleep),]


sum(is.na(sleep$Dream))
# [1] 12

mean(is.na(sleep$Dream))
# [1] 0.1935484

mean(!complete.cases(sleep))
# [1] 0.3225806

# 18.3.1. Tabulating missing values
# install.packages("mice")
library(mice)
data(sleep, package = "VIM")
md.pattern(sleep) # figure 18.2

# 18.3.2. Exploring missing data visually
library("VIM")
aggr(sleep, prop=FALSE, numbers=TRUE) # figure 18.2-2

matrixplot(sleep, sortby = "BrainWgt") # figure 18.3

marginplot(sleep[c("Gest", "Dream")], pch = c(20),
           col = c("darkgray", "red", "blue")) # figure 18.4

# 18.3.3. Using correlations to explore missing values
x <- as.data.frame(abs(is.na(sleep)))

head(sleep, 5)
#    BodyWgt BrainWgt NonD Dream Sleep Span Gest Pred Exp Danger
# 1 6654.000   5712.0   NA    NA   3.3 38.6  645    3   5      3
# 2    1.000      6.6  6.3   2.0   8.3  4.5   42    3   1      3
# 3    3.385     44.5   NA    NA  12.5 14.0   60    1   1      1
# 4    0.920      5.7   NA    NA  16.5   NA   25    5   2      3
# 5 2547.000   4603.0  2.1   1.8   3.9 69.0  624    3   5      4

head(x, 5)
#   BodyWgt BrainWgt NonD Dream Sleep Span Gest Pred Exp Danger
# 1       0        0    1     1     0    0    0    0   0      0
# 2       0        0    0     0     0    0    0    0   0      0
# 3       0        0    1     1     0    0    0    0   0      0
# 4       0        0    1     1     0    1    0    0   0      0
# 5       0        0    0     0     0    0    0    0   0      0

y <- x[which(apply(x, 2, sum) > 0)]

cor(y)
#              NonD       Dream       Sleep        Span        Gest
# NonD   1.00000000  0.90711474  0.48626454  0.01519577 -0.14182716
# Dream  0.90711474  1.00000000  0.20370138  0.03752394 -0.12865350
# Sleep  0.48626454  0.20370138  1.00000000 -0.06896552 -0.06896552
# Span   0.01519577  0.03752394 -0.06896552  1.00000000  0.19827586
# Gest  -0.14182716 -0.12865350 -0.06896552  0.19827586  1.00000000

cor(sleep, y, use = "pairwise.complete.obs")
#                 NonD       Dream        Sleep        Span        Gest
# BodyWgt   0.22682614  0.22259108  0.001684992 -0.05831706 -0.05396818
# BrainWgt  0.17945923  0.16321105  0.007859438 -0.07921370 -0.07332961
# NonD              NA          NA           NA -0.04314514 -0.04553485
# Dream    -0.18895206          NA -0.188952059  0.11699247  0.22774685
# Sleep    -0.08023157 -0.08023157           NA  0.09638044  0.03976464
# Span      0.08336361  0.05981377  0.005238852          NA -0.06527277
# Gest      0.20239201  0.05140232  0.159701523 -0.17495305          NA
# Pred      0.04758438 -0.06834378  0.202462711  0.02313860 -0.20101655
# Exp       0.24546836  0.12740768  0.260772984 -0.19291879 -0.19291879
# Danger    0.06528387 -0.06724755  0.208883617 -0.06666498 -0.20443928


#=======================================================================
# 18.6. Complete-case analysis(listwise deletion)
options(digits = 1)
cor(na.omit(sleep))
#          BodyWgt BrainWgt NonD Dream Sleep  Span  Gest  Pred  Exp Danger
# BodyWgt     1.00     0.96 -0.4 -0.07  -0.3  0.47  0.71  0.10  0.4   0.26
# BrainWgt    0.96     1.00 -0.4 -0.07  -0.3  0.63  0.73 -0.02  0.3   0.15
# NonD       -0.39    -0.39  1.0  0.52   1.0 -0.37 -0.61 -0.35 -0.6  -0.53
# Dream      -0.07    -0.07  0.5  1.00   0.7 -0.27 -0.41 -0.40 -0.5  -0.57
# Sleep      -0.34    -0.34  1.0  0.72   1.0 -0.38 -0.61 -0.40 -0.6  -0.60
# Span        0.47     0.63 -0.4 -0.27  -0.4  1.00  0.65 -0.17  0.3   0.01
# Gest        0.71     0.73 -0.6 -0.41  -0.6  0.65  1.00  0.09  0.6   0.31
# Pred        0.10    -0.02 -0.4 -0.40  -0.4 -0.17  0.09  1.00  0.6   0.93
# Exp         0.41     0.32 -0.6 -0.50  -0.6  0.32  0.57  0.63  1.0   0.79
# Danger      0.26     0.15 -0.5 -0.57  -0.6  0.01  0.31  0.93  0.8   1.00

# alternative
# cor(sleep, use="complete.obs)


#========================================================================
fit <- lm(Dream ~ Span + Gest, data = na.omit(sleep))
summary(fit)
# Call:
# lm(formula = Dream ~ Span + Gest, data = na.omit(sleep))
# 
# Residuals:
#    Min     1Q Median     3Q    Max 
# -2.333 -0.915 -0.221  0.382  4.183 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.480122   0.298476    8.31  3.7e-10 ***
# Span        -0.000472   0.013130   -0.04    0.971    
# Gest        -0.004394   0.002081   -2.11    0.041 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1 on 39 degrees of freedom
# Multiple R-squared:  0.167, Adjusted R-squared:  0.125 
# F-statistic: 3.92 on 2 and 39 DF,  p-value: 0.0282


# 18.7. Multiple imputation
library(mice)
data(sleep, package = "VIM")
imp <- mice(sleep, seed = 1234)

fit <- with(imp, lm(Dream ~ Span + Gest))
pooled <- pool(fit)
summary(pooled)
#          term estimate std.error statistic df p.value
# 1 (Intercept)    2.597     0.249      10.4 52   2e-14
# 2        Span   -0.004     0.012      -0.3 56   7e-01
# 3        Gest   -0.004     0.001      -3.0 55   5e-03



imp
# Class: mids
# Number of multiple imputations:  5 
# Imputation methods:
#   BodyWgt BrainWgt     NonD    Dream    Sleep     Span     Gest     Pred      Exp   Danger 
#        ""       ""    "pmm"    "pmm"    "pmm"    "pmm"    "pmm"       ""       ""       "" 
# PredictorMatrix:
#          BodyWgt BrainWgt NonD Dream Sleep Span Gest Pred Exp Danger
# BodyWgt        0        1    1     1     1    1    1    1   1      1
# BrainWgt       1        0    1     1     1    1    1    1   1      1
# NonD           1        1    0     1     1    1    1    1   1      1
# Dream          1        1    1     0     1    1    1    1   1      1
# Sleep          1        1    1     1     0    1    1    1   1      1
# Span           1        1    1     1     1    0    1    1   1      1
# Number of logged events:  5 
#   it im  dep meth   out
# 1  3  2 Span  pmm Sleep
# 2  3  2 Gest  pmm Sleep
# 3  4  2 Span  pmm Sleep
# 4  4  2 Gest  pmm Sleep
# 5  4  4 Span  pmm Sleep

imp$imp$Dream
#      1   2   3   4   5
# 1  0.0 0.5 0.5 0.5 0.3
# 3  0.5 1.4 1.5 1.5 1.3
# 4  3.6 4.1 3.1 4.1 2.7
# 14 0.3 1.0 0.5 0.0 0.0
# 24 3.6 0.8 1.4 1.4 0.9
# 26 2.4 0.5 3.9 3.4 1.2
# 30 2.6 0.8 2.4 2.2 3.1
# 31 0.6 1.3 1.2 1.8 2.1
# 47 1.3 1.8 1.8 1.8 3.9
# 53 0.5 0.5 0.6 0.5 0.3
# 55 2.6 3.6 2.4 1.8 0.5
# 62 1.5 3.4 3.9 3.4 2.2
str(imp)

options(digits = 3)
dataset3 <- complete(imp, action = 3)
head(dataset3,10)


# 18.8.1. Pairwise deletion
cor(sleep, use = "pairwise.complete.obs")
#          BodyWgt BrainWgt   NonD  Dream  Sleep    Span   Gest    Pred    Exp  Danger
# BodyWgt   1.0000   0.9342 -0.376 -0.109 -0.307  0.3025  0.651  0.0595  0.338  0.1336
# BrainWgt  0.9342   1.0000 -0.369 -0.105 -0.358  0.5093  0.747  0.0339  0.368  0.1459
# NonD     -0.3759  -0.3692  1.000  0.514  0.963 -0.3844 -0.595 -0.3182 -0.544 -0.4839
# Dream    -0.1094  -0.1051  0.514  1.000  0.727 -0.2957 -0.451 -0.4475 -0.537 -0.5793
# Sleep    -0.3072  -0.3581  0.963  0.727  1.000 -0.4102 -0.631 -0.3958 -0.642 -0.5877
# Span      0.3025   0.5093 -0.384 -0.296 -0.410  1.0000  0.615 -0.1025  0.360  0.0618
# Gest      0.6511   0.7472 -0.595 -0.451 -0.631  0.6148  1.000  0.2005  0.638  0.3786
# Pred      0.0595   0.0339 -0.318 -0.447 -0.396 -0.1025  0.201  1.0000  0.618  0.9160
# Exp       0.3383   0.3678 -0.544 -0.537 -0.642  0.3604  0.638  0.6182  1.000  0.7872
# Danger    0.1336   0.1459 -0.484 -0.579 -0.588  0.0618  0.379  0.9160  0.787  1.0000
```
