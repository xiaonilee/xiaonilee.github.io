---
title: "Chapter 12. Resampling statistics and bootstrapping"
date: 2020-11-05
lastmod: 2020-11-05
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

Notebook of Reading Books: R in Action_Chapter 12.

<!--more-->

## This chapter covers

- Understanding the logic of permutation tests

- Applying permutation tests to linear models

- Using bootstrapping to obtain confidence intervals

### 12.1. Permutation tests

- randomization or re-randomization tests

### 12.2. Permutation test with the coin package

#### 12.2.1. Independent two-sample and k-sample tests

- `wilcox.test()`
- `oneway_test()`

#### 12.2.2. Independence in contingency tables

- `chisq_test()` or `cmh_test()`

#### 12.2.3. Independence between numeric variables

- `spearman_test()`

#### 12.2.4. Dependent two-sample and k-sample tests

- `wilcoxsign_test()`

### 12.3. Permutation tests with the lmPerm package

#### 12.3.1. Simple and polynomial regression

- Using `lmp()` instead of `lm()`

#### 12.3.2. Multiple regression

- Applying the `lmp()` function

#### 12.3.3. One-way ANOVA and ANCOVA

- `aovp()`

#### 12.3.4. Two-way ANOVA

- `aovp()`
  
### 12.4. Additional comments on permutation tests

### 12.5. Bootstrapping

### 12.6. Bootstrapping with the boot package

- `install.packages("boot")`

#### 12.6.1. Bootstrapping a single statistic

- Figure 12.1. Distribution of bootstrapped **R-squared** values.

  ![fig121](fig121.png)

#### 12.6.2. Bootstrapping several statistics

- Figure 12.2. Distribution of bootstrapping `regression coefficients` for **car weight**.

  ![fig122](fig122.png)

- Figure 12.3. Distribution of bootstrapping `regression coefficients` for **engine displacement**.

  ![fig123](fig123.png)

Attach is the [Script](chapter12.R) of chapter12.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 12.2.1. Independent two-sample and k-sample tests
# code listing 12.1. t-test versus one-way permutation test for the hypothetical data
library(coin)
score <- c(40, 57, 45, 55, 58, 57, 64, 55, 62, 65)
treatment <- factor(c(rep("A", 5), rep("B", 5)))
mydata <- data.frame(treatment, score)
t.test(score~treatment, data = mydata, var.equal=T)
#         Two Sample t-test
# 
# data:  score by treatment
# t = -2, df = 8, p-value = 0.05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -19.04  -0.16
# sample estimates:
# mean in group A mean in group B 
#               51              61

oneway_test(score~treatment, data = mydata, distribution="exact")
#         Exact Two-Sample Fisher-Pitman Permutation Test
# 
# data:  score by treatment (A, B)
# Z = -2, p-value = 0.07
# alternative hypothesis: true mu is not equal to 0

library(MASS)
UScrime <- transform(UScrime, So = factor(So))
wilcox_test(Prob ~ So, data = UScrime, distribution="exact")
#         Exact Wilcoxon-Mann-Whitney Test
# 
# data:  Prob by So (0, 1)
# Z = -4, p-value = 8e-05
# alternative hypothesis: true mu is not equal to 0


library(multcomp)
set.seed(1234)
oneway_test(response~trt, data = cholesterol,
            distribution=approximate(nresample=9999))

#         Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  response by trt (1time, 2times, 4times, drugD, drugE)
# chi-squared = 36, p-value <1e-04

#===============================================================
# 12.2.2. Independence in contingency tables
library(coin)
library(vcd)
Arthritis <- transform(Arthritis,
              Improved=as.factor(as.numeric(Improved)))
set.seed(1234)
chisq_test(Treatment~Improved, data = Arthritis,
           distribution=approximate(nresample=9999))

#         Approximative Pearson Chi-Squared Test
# 
# data:  Treatment by Improved (1, 2, 3)
# chi-squared = 13, p-value = 0.001


Arthritis <- transform(Arthritis,
                       Improved=as.factor(Improved))
set.seed(1234)
chisq_test(Treatment~Improved, data = Arthritis,
           distribution=approximate(nresample=9999))
#         Approximative Pearson Chi-Squared Test
# 
# data:  Treatment by Improved (1, 2, 3)
# chi-squared = 13, p-value = 0.002


#===============================================
# 12.2.3. Independence between numeric variables

states <- as.data.frame(state.x77)
set.seed(1234)
spearman_test(Illiteracy~Murder, data=states,
               distribution=approximate(nresample=9999))
#         Approximative Spearman Correlation Test
# 
# data:  Illiteracy by Murder
# Z = 5, p-value <1e-04
# alternative hypothesis: true rho is not equal to 0


#====================================================
# 12.2.4. Dependent two-sample and k-sample tests
library(coin)
library(MASS)
wilcoxsign_test(U1~U2, data = UScrime, distribution="exact")

#         Exact Wilcoxon-Pratt Signed-Rank Test
# 
# data:  y by x (pos, neg) 
# stratified by block
# Z = 6, p-value = 1e-14
# alternative hypothesis: true mu is not equal to 0

#====================================================
# 12.3.1. Simple and polynomial regression
# code listing 12.2. Permutation tests for simple linear regression
# install.packages("lmPerm")
library(lmPerm)
set.seed(1234)
fit <- lmp(weight~height, data=women, perm="Prob")
# [1] "Settings:  unique SS : numeric variables centered"

summary(fit)
# Call:
# lmp(formula = weight ~ height, data = women, perm = "Prob")
# 
# Residuals:
#    Min     1Q Median     3Q    Max 
# -1.733 -1.133 -0.383  0.742  3.117 
# 
# Coefficients:
#        Estimate Iter Pr(Prob)    
# height     3.45 5000   <2e-16 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.5 on 13 degrees of freedom
# Multiple R-Squared: 0.991, Adjusted R-squared: 0.99 
# F-statistic: 1.43e+03 on 1 and 13 DF,  p-value: 1.09e-14


#===============================================================
# code listing 12.3. Permutation tests for polynomial regression

library(lmPerm)
set.seed(1234)
fit <- lmp(weight~height + I(height^2), data=women, perm="Prob")
# [1] "Settings:  unique SS : numeric variables centered"
summary(fit)
# Call:
# lmp(formula = weight ~ height + I(height^2), data = women, perm = "Prob")
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.50941 -0.29611 -0.00941  0.28615  0.59706 
# 
# Coefficients:
#             Estimate Iter Pr(Prob)    
# height       -7.3483 5000   <2e-16 ***
#   I(height^2)   0.0831 5000   <2e-16 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.38 on 12 degrees of freedom
# Multiple R-Squared: 0.999, Adjusted R-squared: 0.999 
# F-statistic: 1.14e+04 on 2 and 12 DF,  p-value: <2e-16


#=======================================================
# 12.3.2. Multiple regression
# code listing 12.4. Permutation tests for multiple regression
library(lmPerm)
set.seed(1234)
states <- as.data.frame(state.x77)
fit <- lmp(Murder~Population + Illiteracy+Income+Frost,
           data=states, perm="Prob")
# [1] "Settings:  unique SS : numeric variables centered"

summary(fit)
# Call:
# lmp(formula = Murder ~ Population + Illiteracy + Income + Frost, 
#       data = states, perm = "Prob")
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -4.7960 -1.6495 -0.0811  1.4815  7.6210 
# 
# Coefficients:
#            Estimate Iter Pr(Prob)    
# Population 2.24e-04   51   1.0000    
# Illiteracy 4.14e+00 5000   0.0004 ***
# Income     6.44e-05   51   1.0000    
# Frost      5.81e-04   51   0.8627    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.5 on 45 degrees of freedom
# Multiple R-Squared: 0.567, Adjusted R-squared: 0.528 
# F-statistic: 14.7 on 4 and 45 DF,  p-value: 9.13e-08

#==============================================================
# 12.3.3. One-way ANOVA and ANCOVA
# code listing 12.5. Permutation test for One-Way ANOVA
library(lmPerm)
library(multcomp)
set.seed(1234)
fit <- aovp(response~trt, data=cholesterol, perm="Prob")
# [1] "Settings:  unique SS "
summary(fit)
# Component 1 :
#             Df R Sum Sq R Mean Sq Iter  Pr(Prob)    
# trt          4  1351.37    337.84 5000 < 2.2e-16 ***
# Residuals   45   468.75     10.42                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#==============================================================
# code listing 12.6. Permutation test for one-way ANCOVA
library(lmPerm)
set.seed(1234)
fit <- aovp(weight ~ gesttime + dose, data=litter, perm="Prob")
# [1] "Settings:  unique SS : numeric variables centered"

summary(fit)
# Component 1 :
#             Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# gesttime     1   161.49   161.493 5000   0.0006 ***
# dose         3   137.12    45.708 5000   0.0392 *  
# Residuals   69  1151.27    16.685                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#==============================================================
# 12.3.4. Two-way ANOVA
# code listing 12.7. Permutation test for two-way ANOVA
library(lmPerm)
set.seed(1234)
fit <- aovp(len~supp*dose, data=ToothGrowth, perm="Prob")
# [1] "Settings:  unique SS : numeric variables centered"
summary(fit)
# Component 1 :
#             Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# supp         1   205.35    205.35 5000  < 2e-16 ***
# dose         1  2224.30   2224.30 5000  < 2e-16 ***
# supp:dose    1    88.92     88.92 2032  0.04724 *  
# Residuals   56   933.63     16.67                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#==============================================================
# 12.6.1. Bootstrapping a single statistic
rsq <- function(formula, data, indices) {
  d <- data[indices,]
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
  
}

library(boot)
set.seed(1234)
results <- boot(data=mtcars, statistic=rsq,
                R=1000, formula=mpg~wt+disp)

print(results)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
# boot(data = mtcars, statistic = rsq, R = 1000, formula = mpg ~ 
#          wt + disp)
# 
# 
# Bootstrap Statistics :
#      original     bias    std. error
# t1* 0.7809306 0.01379126  0.05113904

plot(results) # figure 12.1
boot.ci(results, type = c("perc", "bca"))
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# 
# CALL : 
# boot.ci(boot.out = results, type = c("perc", "bca"))
# 
# Intervals : 
# Level     Percentile            BCa          
# 95%   ( 0.6753,  0.8835 )   ( 0.6344,  0.8561 )  
# Calculations and Intervals on Original Scale
# Some BCa intervals may be unstable


#=============================================
# 12.6.2. Bootstrapping several statistics
bs <- function(formula, data, indices) {
  d <- data[indices,]
  fit <- lm(formula, data=d)
  return(coef(fit))
}

library(boot)
set.seed(1234)
results <- boot(data=mtcars, statistic=bs,
                R=1000, formula=mpg~wt+disp)
print(results)
# ORDINARY NONPARAMETRIC BOOTSTRAP
# 
# 
# Call:
# boot(data = mtcars, statistic = bs, R = 1000, formula = mpg ~ 
#          wt + disp)
# 
# 
# Bootstrap Statistics :
#        original        bias    std. error
# t1* 34.96055404  4.715497e-02 2.546106756
# t2* -3.35082533 -4.908125e-02 1.154800744
# t3* -0.01772474  6.230927e-05 0.008518022

plot(results, index = 2)  # figure 12.2

boot.ci(results, type = "bca", index = 2) # t2: car weight
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# 
# CALL : 
# boot.ci(boot.out = results, type = "bca", index = 2)
# 
# Intervals : 
#   Level       BCa          
# 95%   (-5.477, -0.937 )  
# Calculations and Intervals on Original Scale

boot.ci(results, type = "bca", index = 3) # t3: engine displacement
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# 
# CALL : 
# boot.ci(boot.out = results, type = "bca", index = 3)
# 
# Intervals : 
#   Level       BCa          
# 95%   (-0.0334, -0.0011 )  
# Calculations and Intervals on Original Scale

plot(results, index = 3)  # figure 12.3
```
