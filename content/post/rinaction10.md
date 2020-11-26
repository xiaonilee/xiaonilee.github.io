---
title: "Chapter 10. Power analysis"
date: 2020-11-02
lastmod: 2020-11-02
draft: false
tags: ["R", "R in Action", "Bioinformatics", "Book"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

Notebook of Reading Books: R in Action_Chapter 10.

<!--more-->

## This chapter covers

- Determining sample size requirements

- Calculating effect sizes

- Assessing statistical power

### 10.1. A quick review of hypothesis testing

$$ H_{0} : u_{1} - u_{2} = 0 $$ 
$$ H_{1} : u_{1} - u_{2} ≠ 0 $$

- Based on the sample data, you can calculate the statistic.
  
  ![tab1](tab1.png)

- Figure 10.1. Four primary quantities considered in a study design power analysis. 

- Given any three, you can calculate the fourth.

  ![fig101](fig101.png)

### 10.2. Implementing power analysis with the pwr package

  ![tab101](tab101.png)

#### 10.2.1. t-tests

- equal n : pwr.t.test(n=, d=, sig.level=, power=, alternative=)

- unequal n : pwr.t2n.test(n1=, n2=, d=, sig.level=, power=, alternative=)

#### 10.2.2. ANOVA

- pwr.anova.test(k=, n=, f=, sig.level=, power=)

#### 10.2.3. Correlations

- pwr.r.test(n=, r=, sig.level=, power=, alternative=)

#### 10.2.4. Linear models

- pwr.f2.test(u=, v=, f2=, sig.level=, power=)

#### 10.2.5. Tests of proportions

- equal n: pwr.2p.test(h=, n=, sig.level=, power=)

- unequal n: pwr.2p2n.test(h =, n1 =, n2 =, sig.level=, power=)

#### 10.2.6. Chi-square tests

- pwr.chisq.test(w =, N = , df = , sig.level =, power =)

#### 10.2.7. Choosing an appropriate effect size in novel situations

- Without experiences for the reference
  
  ![tab103](tab103.png)

- Figure 10.2. Sample size needed to detect various effect sizes in a one-way ANOVA with five groups (assuming a power of 0.90 and significance level of 0.05).

  ![fig102](fig102.png)

### 10.3. Creating power analysis plots

- Figure 10.3. Sample size curves for detecting a significant correlation at various power levels.

  ![fig103](fig103.png)

### 10.4. Other packages

- The last five are particularly focused on power analysis in **genetic studies**. 

- Genome-wide association studies (`GWAS`) are studies used to identify genetic associations with observable traits.

  ![tab104](tab104.png)

Attach is the [Script](chapter10.R) of chapter10.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 10.2.1. t-tests
# Question 1
# install.packages("pwr")
library(pwr)
pwr.t.test(d=0.8, sig.level = 0.05, power=0.9, type = "two.sample",
           alternative = "two.sided")
#     Two-sample t test power calculation 
# 
#              n = 33.82555
#              d = 0.8
#       sig.level = 0.05
#           power = 0.9
#     alternative = two.sided
# 
# NOTE: n is number in *each* group


#===================================================================
# question 2
pwr.t.test(n=20, d=0.5, sig.level = 0.01, type = "two.sample",
           alternative = "two.sided")
#     Two-sample t test power calculation 
# 
#              n = 20
#              d = 0.5
#      sig.level = 0.01
#          power = 0.1439551
#    alternative = two.sided
# 
# NOTE: n is number in *each* group

#===================================================================
# 10.2.2. ANOVA
pwr.anova.test(k=5, f=0.25,sig.level = 0.05, power = 0.8)

#     Balanced one-way analysis of variance power calculation 
# 
#              k = 5
#              n = 39.1534
#              f = 0.25
#      sig.level = 0.05
#          power = 0.8
# 
# NOTE: n is number in each group

#===================================================================
# 10.2.3. Correlations

# question 3
pwr.r.test(r=0.25, sig.level = 0.05, power = 0.90, alternative = "greater")
#     approximate correlation power calculation (arctangh transformation) 
# 
#              n = 133.2803
#              r = 0.25
#      sig.level = 0.05
#          power = 0.9
# alternative = greater

#===================================================================
# 10.2.4. Linear models
pwr.f2.test(u=3, f2=0.0769, sig.level = 0.05, power=0.90)

#       Multiple regression power calculation 
# 
#                u = 3
#                v = 184.2426
#                f2 = 0.00769
#         sig.level = 0.05
#             power = 0.9

#       v=N-K-1-->N=v+K+1=185+7+1=193

#===================================================================
# 10.2.5. Tests of proportions
# question 4

pwr.2p.test(h=ES.h(0.65, 0.6), sig.level = 0.05, power = 0.9,
            alternative = "greater")
#   Difference of proportion power calculation for binomial 
#   distribution (arcsine transformation) 
# 
#            h = 0.1033347
#            n = 1604.007
#    sig.level = 0.05
#        power = 0.9
#  alternative = greater
# 
# NOTE: same sample sizes

#===================================================================
# 10.2.6. Chi-square tests

prob <- matrix(c(.42, .28, .03, .07, .10, .10), byrow=TRUE, nrow=3)
ES.w2(prob)
# [1] 0.1853198
pwr.chisq.test(w=0.1853, df=2, sig.level = 0.05, power=0.90)
#   Chi squared power calculation 
# 
#             w = 0.1853
#             N = 368.5317
#            df = 2
#     sig.level = 0.05
#         power = 0.9
# 
# NOTE: N is the number of observations

#===================================================================
# 10.2.7. Choosing an appropriate effect size in novel situations
# code listing 10.1. Sample sizes for detecting significant effects in a one-way ANOVA
library(pwr)
es <- seq(0.1, 0.5, 0.01)
nes <- length(es)

samsize <- NULL
for (i in 1:nes){
  result <- pwr.anova.test(k=5, f=es[i], sig.level = 0.05, power = 0.9)
  samsize[i] <- ceiling(result$n)
}

plot(samsize, es, type = "l", lwd=2, col="red",
     ylab = "Effect Size",
     xlab = "Sample Size (per cell)",
     main = "One Way ANOVA with Power=0.90, and Alpha=0.05")



# 10.3. Creating power analysis plots
# code listing Sample size curves for detecting correlations of various sizes

library(pwr)
r <- seq(0.1, 0.5, 0.01)
nr <- length(r)
nr

p <- seq(0.4, 0.9, 0.1)
np <- length(p)
np

samsize <- array(numeric(nr*np), dim = c(nr, np))
for (i in 1:np){
  for (j in 1:nr){
    result <- pwr.r.test(n =NULL, r = r[j],
    sig.level = 0.05, power = p[i],
    alternative = "two.sided")
    samsize[j, i] <- ceiling(result$n)
  }
}

xrange <- range(r)
yrange <- round(range(samsize))
colors <- rainbow(length(p))
plot(xrange, yrange, type = "n",
     xlab = "Correlation Coefficient (r)",
     ylab = "Sample Size(n")

for (i in 1:np){
  lines(r, samsize[,i], type = "l", lwd = 2, col=colors[i])
}

abline(v=0, h=seq(0, yrange[2],50), lty=2, col="grey89")
abline(h=0, v=seq(xrange[2], 0.02), lty=2, col="grey89")
title("Sample Size Estimation for Correlation Studies\n
      Sig = 0.05 (Two-tailed)")
legend("topright", title = "Power", as.character(p),
       fill = colors)
```
