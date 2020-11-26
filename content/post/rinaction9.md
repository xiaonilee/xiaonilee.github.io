---
title: "Chapter 9. Analysis of variance"
date: 2020-10-30
lastmod: 2020-10-30
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

Notebook of Reading Books: R in Action_Chapter 9.

<!--more-->

## This chapter covers

- Using R to model basic experimental designs

- Fitting and interpreting ANOVA type models

- Evaluating model assumptions

### 9.1. A crash course on terminology

- ANOVA
  - F test 

### 9.2. Fitting ANOVA models

#### 9.2.1. The aov() function

![tab94](tab94.png)

#### 9.2.2. The order of formula terms

- with `Type I` (sequential), by default.

- `help(Anova, package="car")`

### 9.3. One-way ANOVA

- Figure 9.1. **Treatment group means** with **95 percent confidence intervals** for five cholesterol-reducing drug regiments.

  ![fig91](fig91.png)

#### 9.3.1. Multiple comparisons

- Figure 9.2. Plot of Tukey HSD **pairwise mean** comparisons.

  ![fig92](fig92.png)

- Figure 9.3. Tukey HSD tests provided by the multcomp package.

  ![fig93](fig93.png)

- Figure 9.4. Test of normality.

  ![fig94](fig94.png)

### 9.4. One-way ANCOVA

#### 9.4.1. Assessing test assumptions

- Assume homogeneity of regression slopes

#### 9.4.2. Visualizing the results

- Figure 9.5. Plot of the relationship between gestation time and birth weight for each of four drug treatment groups.

  ![fig95](fig95.png)

### 9.5. Two-way factorial ANOVA

- Figure 9.6. Interaction between dose and delivery mechanism on tooth growth. The plot of means was created using the `interaction.plot()` function.
  
  ![fig96](fig96.png)

- Figure 9.7. Interaction between dose and delivery mechanism on tooth growth. The mean plot with 95 percent confidence intervals was created by the `plotmeans()` function.

  ![fig97](fig97.png)

- Figure 9.8. Main effects and two-way interaction for the ToothGrowth dataset. This plot was created by the `interaction2way()` function in the "HH" package.

  ![fig98](fig98.png)

### 9.6. Repeated measures ANOVA

- Figure 9.9. Interaction of ambient CO2 concentration and plant type on CO2 uptake. Graph produced by the `interaction.plot()` function.

  ![fig99](fig99.png)

- Figure 9.10. Interaction of ambient CO2 concentration and plant type on CO2 uptake. Graph produced by the `boxplot()` function.

  ![fig910](fig910.png)

### 9.7. Multivariate analysis of variance (MANOVA)

#### 9.7.1. Assessing test assumptions

- Figure 9.11. A Q-Q plot for assessing multivariate normality.

  ![fig911](fig911.png)

- `aq.plot()` in the **mvoutlier** package to check `multivariate outliers`

  ![outliers](outliers.png)
  
#### 9.7.2. Robust MANOVA

- one-way MANOVA
  
  - with `Wilks.test()` in the `rrcov` package.

- a nonparametric MANOVA
  
  - with `adonis()` in the `vegan` package

### 9.8. ANOVA as regression

- Table 9.6. Built-in contrasts
  
  ![Tab96](tab96.png)

Attach is the [Script](chapter9.R) of chapter9.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 9.3. One-way ANOVA
# code listing 9.1
library(multcomp)
attach(cholesterol)
colnames(cholesterol)
table(trt)  # group sample sizes
# trt
#  1time 2times 4times  drugD  drugE 
#     10     10     10     10     10 

aggregate(response, by=list(trt), mean) # group means
#   Group.1        x
# 1   1time  5.78197
# 2  2times  9.22497
# 3  4times 12.37478
# 4   drugD 15.36117
# 5   drugE 20.94752

aggregate(response, by=list(trt), sd) # group standard deviations
#   Group.1        x
# 1   1time 2.878113
# 2  2times 3.483054
# 3  4times 2.923119
# 4   drugD 3.454636
# 5   drugE 3.345003

fit <- aov(response ~ trt)
summary(fit)  # Test for group differences(ANOVA)

#             Df Sum Sq Mean Sq F value   Pr(>F)    
# trt          4 1351.4   337.8   32.43 9.82e-13 ***
# Residuals   45  468.8    10.4                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(gplots)
plotmeans(response ~ trt,
          xlab = "Treatment",
          ylab = "Response",
          main = "Mean Plot \nwith 95% CI") # figure 9.1
detach(cholesterol)

# 9.3.1. Multiple comparisons
# code listing 9.2. Tukey HSD pairwise group comparisons
TukeyHSD(fit)
# Tukey multiple comparisons of means
#   95% family-wise confidence level
# 
# Fit: aov(formula = response ~ trt)
# 
# $trt
#                   diff        lwr       upr     p adj
# 2times-1time   3.44300 -0.6582817  7.544282 0.1380949
# 4times-1time   6.59281  2.4915283 10.694092 0.0003542
# drugD-1time    9.57920  5.4779183 13.680482 0.0000003
# drugE-1time   15.16555 11.0642683 19.266832 0.0000000
# 4times-2times  3.14981 -0.9514717  7.251092 0.2050382
# drugD-2times   6.13620  2.0349183 10.237482 0.0009611
# drugE-2times  11.72255  7.6212683 15.823832 0.0000000
# drugD-4times   2.98639 -1.1148917  7.087672 0.2512446
# drugE-4times   8.57274  4.4714583 12.674022 0.0000037
# drugE-drugD    5.58635  1.4850683  9.687632 0.0030633


par(las=2)
par(mar=c(5,8,4,2))
plot(TukeyHSD(fit)) # figure 9.2


library(multcomp)
par(mar=c(5,4,6,2))
tuk <- glht(fit, linfct=mcp(trt="Tukey"))
plot(cld(tuk, level=.05),col="lightgrey") # figure 9.3


# 9.3.2. Assessing test assumptions
# figure 9.4
library(car)
qqPlot(lm(response ~ trt, data=cholesterol),
         simulate=TRUE, main="Q-Q Plot", labels=FALSE)


# code
bartlett.test(response ~ trt, data=cholesterol)

#         Bartlett test of homogeneity of variances
# 
# data:  response by trt
# Bartlett's K-squared = 0.57975, df = 4, p-value = 0.9653


outlierTest(fit)
# No Studentized residuals with Bonferroni p < 0.05
# Largest |rstudent|:
#   rstudent unadjusted p-value Bonferroni p
# 19 2.251149           0.029422           NA

#=============================================================

# 9.4. One-way ANCOVA
# code listing 9.3. One-way ANCOVA
data("litter", package = "multcomp")
attach(litter)
table(dose)
# dose
#   0   5  50 500 
#   20  19  18  17 

aggregate(weight, by=list(dose), mean)
#   Group.1        x
# 1       0 32.30850
# 2       5 29.30842
# 3      50 29.86611
# 4     500 29.64647

fit <- aov(weight ~ gesttime + dose) 
summary(fit)

#             Df Sum Sq Mean Sq F value   Pr(>F)    
# trt          4 1351.4   337.8   32.43 9.82e-13 ***
# Residuals   45  468.8    10.4                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(effects)
effect("dose", fit)

#  dose effect
# dose
#        0        5       50      500 
# 32.35367 28.87672 30.56614 29.33460 


# code listing 9.4 Multiple comparisons employing user-supplied contrasts
library(multcomp)
contrast <- rbind("no drug vs. drug" = c(3, -1, -1, -1))
summary(glht(fit, linfct=mcp(dose=contrast)))
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: aov(formula = weight ~ gesttime + dose)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)  
# no drug vs. drug == 0    8.284      3.209   2.581    0.012 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

help("glht")

# 9.4.1. Assessing test assumptions
# code listing 9.5. Testing for homogeneity of regression slopes
library(multcomp)
fit2 <- aov(weight ~ gesttime*dose, data = litter)
summary(fit2)

#               Df Sum Sq Mean Sq F value  Pr(>F)   
# gesttime       1  134.3  134.30   8.289 0.00537 **
# dose           3  137.1   45.71   2.821 0.04556 * 
# gesttime:dose  3   81.9   27.29   1.684 0.17889   
# Residuals     66 1069.4   16.20                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#==================================================================
# 9.4.2. Visualizing the results
# install.packages("HH")
library(HH)
ancova(weight ~ gesttime + dose, data=litter) # figure 9.5

# 9.5. Two-way factorial ANOVA
# code listing 9.6 Two-way ANOVA
attach(ToothGrowth)
head(ToothGrowth)
#    len supp dose
# 1  4.2   VC  0.5
# 2 11.5   VC  0.5
# 3  7.3   VC  0.5
# 4  5.8   VC  0.5
# 5  6.4   VC  0.5
# 6 10.0   VC  0.5

class(dose)
colnames(ToothGrowth)
# [1] "len"  "supp" "dose"

table(supp,dose)
#     dose
# supp 0.5  1  2
# OJ  10 10 10
# VC  10 10 10

aggregate(len, by=list(supp, dose), mean)
#   Group.1 Group.2     x
# 1      OJ     0.5 13.23
# 2      VC     0.5  7.98
# 3      OJ     1.0 22.70
# 4      VC     1.0 16.77
# 5      OJ     2.0 26.06
# 6      VC     2.0 26.14

aggregate(len, by=list(supp, dose), sd)
#   Group.1 Group.2        x
# 1      OJ     0.5 4.459709
# 2      VC     0.5 2.746634
# 3      OJ     1.0 3.910953
# 4      VC     1.0 2.515309
# 5      OJ     2.0 2.655058
# 6      VC     2.0 4.797731

fit <- aov(len ~ supp*dose)
summary(fit)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# supp         1  205.4   205.4  12.317 0.000894 ***
# dose         1 2224.3  2224.3 133.415  < 2e-16 ***
# supp:dose    1   88.9    88.9   5.333 0.024631 *  
# Residuals   56  933.6    16.7                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dose <- factor(dose)
fit2 <- aov(len ~ supp*dose)
summary(fit2)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# supp         1  205.4   205.4  15.572 0.000231 ***
# dose         2 2426.4  1213.2  92.000  < 2e-16 ***
# supp:dose    2  108.3    54.2   4.107 0.021860 *  
# Residuals   54  712.1    13.2                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

detach(ToothGrowth)

# figure 9.6
interaction.plot(dose, supp, len, type="b",
                 col=c("red","blue"), pch=c(16, 18),
                 main = "Interaction between Dose and Supplement Type")

# figure 9.7
library(gplots)
plotmeans(len ~ interaction(supp, dose, sep=" "),
          connect=list(c(1,3,5),c(2,4,6)),
          col=c("red", "darkgreen"),
          main = "Interaction Plot with 95% CIs",
          xlab="Treatment and Dose Combination")

# figure 9.8
library(HH)
interaction2wt(len~supp*dose)

# 9.6. Repeated measures ANOVA
# code listing 9.7. Repeated measures ANOVA with one between- and within-groups factor
CO2$conc <- factor(CO2$conc)
w1b1 <- subset(CO2, Treatment=='chilled')
fit <- aov(uptake ~ conc*Type + Error(Plant/(conc)), w1b1)
summary(fit)
# Error: Plant
#           Df Sum Sq Mean Sq F value  Pr(>F)   
# Type       1 2667.2  2667.2   60.41 0.00148 **
# Residuals  4  176.6    44.1                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: Plant:conc
#           Df Sum Sq Mean Sq F value   Pr(>F)    
# conc       6 1472.4  245.40   52.52 1.26e-12 ***
# conc:Type  6  428.8   71.47   15.30 3.75e-07 ***
# Residuals 24  112.1    4.67                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

par(las=2)
par(mar=c(10, 4, 4, 2))
# figure 9.9
with(w1b1, interaction.plot(conc,Type,uptake,
                            type="b", col=c("red","blue"), pch=c(16,18),
                            main="Interaction Plot for Plant Type and Concentration"))

# figure 9.10
boxplot(uptake ~ Type*conc, data=w1b1, col=(c("gold", "green")),
        main="Chilled Quebec and Mississippi Plants",
        ylab="Carbon dioxide uptake rate (umol/m^2 sec)")


# 9.7. Multivariate analysis of variance (MANOVA)
# code listing 9.8. One-way MANOVA
library(MASS)
attach(UScereal)
shelf <- factor(shelf)
y <- cbind(calories, fat, sugars)
aggregate(y, by=list(shelf), mean)
#   Group.1 calories       fat    sugars
# 1       1 119.4774 0.6621338  6.295493
# 2       2 129.8162 1.3413488 12.507670
# 3       3 180.1466 1.9449071 10.856821

cov(y)
#            calories       fat     sugars
# calories 3895.24210 60.674383 180.380317
# fat        60.67438  2.713399   3.995474
# sugars    180.38032  3.995474  34.050018

fit <- manova(y ~ shelf)
summary(fit)  
#           Df Pillai approx F num Df den Df    Pr(>F)    
# shelf      2 0.4021   5.1167      6    122 0.0001015 ***
# Residuals 62                                            
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
summary.aov(fit)
# Response calories :
#             Df Sum Sq Mean Sq F value    Pr(>F)    
# shelf        2  50435 25217.6  7.8623 0.0009054 ***
# Residuals   62 198860  3207.4                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response fat :
#             Df Sum Sq Mean Sq F value  Pr(>F)  
# shelf        2  18.44  9.2199  3.6828 0.03081 *
# Residuals   62 155.22  2.5035                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response sugars :
#             Df  Sum Sq Mean Sq F value   Pr(>F)   
# shelf        2  381.33 190.667  6.5752 0.002572 **
# Residuals   62 1797.87  28.998                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# 9.7.1. Assessing test assumptions
# code listing 9.9. Assessing multivariate normality
center <- colMeans(y)
n <- nrow(y)
p <- ncol(y)
cov <- cov(y)
d <- mahalanobis(y, center, cov)

# figure 9.11
coord <- qqplot(qchisq(ppoints(n), df=p),
                d, main="Q-Q Plot Assessing Multivariate Normality",
                ylab="Mahalanobis D2")
abline(a=0,b=1)
identify(coord$x, coord$y, labels = row.names(UScereal))

#================================================================
# install.packages("mvoutlier")
library(mvoutlier)
outliers <- aq.plot(y)
outliers

# 9.7.2. Robust MANOVA
# code listing 9.10. Robust one-way MANOVA
library(rrcov)
Wilks.test(y, shelf, method = "mcd")
#         Robust One-way MANOVA (Bartlett Chi2)
# 
# data:  x
# Wilks' Lambda = 0.51073, Chi2-Value = 25.3745, DF = 5.1632, p-value = 0.0001377
# sample estimates:
#   calories       fat    sugars
# 1 119.8210 0.7010828  5.663143
# 2 128.0407 1.1849576 12.537533
# 3 160.8604 1.6524559 10.352646


# 9.8. ANOVA as regression
library(multcomp)
levels(cholesterol$trt)
# [1] "1time"  "2times" "4times" "drugD"  "drugE"

fit.aov <- aov(response ~ trt, data=cholesterol)
summary(fit.aov)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# trt          4 1351.4   337.8   32.43 9.82e-13 ***
# Residuals   45  468.8    10.4                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# code listing 9.11. A regression approach to the ANOVA problem in section 9.3
fit.lm <- lm(response ~ trt, data = cholesterol)
summary(fit.lm)

# Call:
#   lm(formula = response ~ trt, data = cholesterol)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6.5418 -1.9672 -0.0016  1.8901  6.6008 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    5.782      1.021   5.665 9.78e-07 ***
#   trt2times      3.443      1.443   2.385   0.0213 *  
#   trt4times      6.593      1.443   4.568 3.82e-05 ***
#   trtdrugD       9.579      1.443   6.637 3.53e-08 ***
#   trtdrugE      15.166      1.443  10.507 1.08e-13 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.227 on 45 degrees of freedom
# Multiple R-squared:  0.7425, Adjusted R-squared:  0.7196 
# F-statistic: 32.43 on 4 and 45 DF,  p-value: 9.819e-13

contrasts(cholesterol$trt)
#        2times 4times drugD drugE
# 1time       0      0     0     0
# 2times      1      0     0     0
# 4times      0      1     0     0
# drugD       0      0     1     0
# drugE       0      0     0     1
```
