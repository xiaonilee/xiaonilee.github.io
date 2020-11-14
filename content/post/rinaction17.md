---
title: "Chapter 17. Classification"
date: 2020-11-13
lastmod: 2020-11-13
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

Notebook of Reading Books: R in Action_Chapter 17.

<!--more-->

## This chapter covers

- Classifying with decision trees

- Ensemble classification with random forests

- Creating a support vector machine

- Evaluating classification accuracy

### 17.1 Preparing the data

- Data from the UCI Machine Learning Server
  - <http://archive.ics.uci.edu/ml/index.php>
  - Breast Cancer Wisconsin (Original) Data Set
    - <http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+%28Original%29>

### 17.2 Logistic regression

- `glm()` & `predict()`

### 17.3. Decision trees

#### 17.3.1. Classical decision trees

- Figure 17.1. Complexity parameter vs. cross-validated error. 
  - The dotted line is the upper limit of the one standard deviation rule `(0.12 + 1 * 0.026 = .15)`
  - The plot suggests selecting the tree with the leftmost cp value below the line.

  ![fig171](fig171.png)

- Figure 17.2. Traditional (pruned) decision tree for predicting cancer status. 
  - Start at the top of the tree, moving left if a condition is true or right otherwise. 
  - When an observation hits a terminal node, it’s classified. Each node contains the probability of the classes in that node, along with the percentage of the sample.
  
  ![fig172](fig172.png)

#### 17.3.2 Conditional inference trees

- Figure 17.3. Conditional inference tree for the breast cancer data.

  ![fig173](fig173.png)

### 17.4 Random forests

### 17.5. Support vector machines

#### 17.6. Choosing a best predictive solution

  ![tab171](tab171.png)

- define performance() function

### 17.7 Using the rattle package for data mining

- Installing rattle on MacOS 10.11 (or above) with [R Script](install_rattle.R).

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Installing rattle on MacOS 10.11 (or above)

# install gtk+ 2.24.32_3
system('brew install gtk+')


# macOS version:catalina 
local({
  if (Sys.info()[['sysname']] != 'Darwin') return()
  
  .Platform$pkgType = 'mac.binary.catalina'
  unlockBinding('.Platform', baseenv())
  assign('.Platform', .Platform, 'package:base')
  lockBinding('.Platform', baseenv())
  
  options(
    pkgType = 'both', install.packages.compile.from.source = 'always',
    repos = 'https://macos.rbind.io'
  )
})

install.packages(c('RGtk2', 'cairoDevice', 'rattle'))
```

  ![rattleUI](rattleUI.png)

- Load and setup data: 

  ![rattleUI2](rattleUI20.png)

- Execute

  ![rattleUI22](rattleUI2.png)

- model(up) and draw(bottom)

  ![model](model4.png)
  ![fig179](fig179.png)

- Evaluate
  
  ![evaluate5](evaluate5.png)

Attach is the [Script](chapter17.R) of chapter17.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# Prerequisites
pkgs <- c("rpart", "rpart.plot", "party",
          "randomForest", "e1071")
install.packages(pkgs, depend=TRUE)


# 17.1 Preparing the data
# code listing 17.1. Preparing the breast cancer data
loc <- "http://archive.ics.uci.edu/ml/machine-learning-databases/"
ds <- "breast-cancer-wisconsin/breast-cancer-wisconsin.data"
url <- paste(loc, ds, sep = "")
url


breast <- read.table(url, sep = ",", header = FALSE, na.strings = "?")
head(breast)
#        V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11
# 1 1000025  5  1  1  1  2  1  3  1   1   2
# 2 1002945  5  4  4  5  7 10  3  2   1   2
# 3 1015425  3  1  1  1  2  2  3  1   1   2
# 4 1016277  6  8  8  1  3  4  3  7   1   2
# 5 1017023  4  1  1  3  2  1  3  1   1   2
# 6 1017122  8 10 10  8  7 10  9  7   1   4


names(breast) <- c("ID", "clumpThickness", "sizeUniformity",
                   "shapeUniformity", "maginalAdhesion",
                   "singleEpithelialCellSize", "bareNuclei",
                   "blandChromatin", "normalNucleoli", "mitosis", "class")

table(breast$class)
#   2   4 
# 458 241

df <- breast[-1]

nrow(df)
# [1] 699

df$class <- factor(df$class, levels = c(2, 4),
                   labels = c("benign", "malignant"))
set.seed(1234)
train <- sample(nrow(df), 0.7*nrow(df))
df.train <- df[train,]
df.validate <- df[-train,]
table(df.train$class)
# benign malignant 
#   319       170

table(df.validate$class)
# benign malignant 
#   139        71


#=========================================================
# 17.2 Logistic regression
# code listing 17.2. Logistic regression with glm()
## Fits the logistic regression with df.train dataset
fit.logit <- glm(class ~ ., data = df.train, family = binomial())

##Examines the model
summary(fit.logit)
# Call:
# glm(formula = class ~ ., family = binomial(), data = df.train)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -2.24605  -0.08012  -0.03110   0.00266   2.11576  
# 
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -12.4412     2.0547  -6.055  1.4e-09 ***
# clumpThickness             0.7407     0.2262   3.275  0.00106 ** 
# sizeUniformity            -0.0320     0.3399  -0.094  0.92500    
# shapeUniformity            0.2073     0.3715   0.558  0.57680    
# maginalAdhesion            0.5194     0.1708   3.041  0.00236 ** 
# singleEpithelialCellSize  -0.3217     0.2613  -1.231  0.21831    
# bareNuclei                 0.5851     0.1881   3.111  0.00187 ** 
# blandChromatin             0.8599     0.2923   2.942  0.00326 ** 
# normalNucleoli             0.4036     0.1828   2.208  0.02725 *  
# mitosis                    0.8923     0.3552   2.512  0.01200 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 621.04  on 477  degrees of freedom
# Residual deviance:  52.39  on 468  degrees of freedom
#   (11 observations deleted due to missingness)
# AIC: 72.39
# 
# Number of Fisher Scoring iterations: 9

## Classifies new cases with df.validate dataset
prob <- predict(fit.logit, df.validate, type = "response")
logit.pred <- factor(prob > .5, levels = c(FALSE, TRUE),
                     labels = c("benign", "malignant"))
summary(logit.pred)
#    benign malignant      NA's 
#       130        75         5 


## Evaluates the predictive accuracy
logit.perf <- table(df.validate$class, logit.pred,
                    dnn = c("Actual", "Predicted"))
logit.perf
#            Predicted
# Actual      benign malignant
#   benign       129         6
#   malignant      1        69

summary(logit.perf)
# Number of cases in table: 205 
# Number of factors: 2 
# Test for independence of all factors:
#   Chisq = 176.04, df = 1, p-value = 3.55e-40
(129+69)/(129+6+1+69)*100
# [1] 96.58537



logit.fit.reduced <- step(fit.logit)

# Start:  AIC=72.39
# class ~ clumpThickness + sizeUniformity + shapeUniformity + maginalAdhesion + 
#   singleEpithelialCellSize + bareNuclei + blandChromatin + 
#   normalNucleoli + mitosis
# 
# Df Deviance    AIC
# - sizeUniformity            1   52.399 70.399
# - shapeUniformity           1   52.686 70.686
# - singleEpithelialCellSize  1   53.910 71.910
# <none>                          52.390 72.390
# - normalNucleoli            1   57.465 75.465
# - mitosis                   1   57.966 75.966
# - blandChromatin            1   62.856 80.856
# - maginalAdhesion           1   63.044 81.044
# - bareNuclei                1   64.982 82.982
# - clumpThickness            1   68.323 86.323
# 
# Step:  AIC=70.4
# class ~ clumpThickness + shapeUniformity + maginalAdhesion + 
#   singleEpithelialCellSize + bareNuclei + blandChromatin + 
#   normalNucleoli + mitosis
# 
# Df Deviance    AIC
# - shapeUniformity           1   52.876 68.876
# - singleEpithelialCellSize  1   53.918 69.918
# <none>                          52.399 70.399
# - normalNucleoli            1   57.916 73.916
# - mitosis                   1   58.024 74.024
# - blandChromatin            1   63.272 79.272
# - maginalAdhesion           1   63.735 79.735
# - bareNuclei                1   64.985 80.985
# - clumpThickness            1   68.728 84.728
# 
# Step:  AIC=68.88
# class ~ clumpThickness + maginalAdhesion + singleEpithelialCellSize + 
#   bareNuclei + blandChromatin + normalNucleoli + mitosis
# 
# Df Deviance    AIC
# - singleEpithelialCellSize  1   54.154 68.154
# <none>                          52.876 68.876
# - mitosis                   1   59.402 73.402
# - normalNucleoli            1   60.929 74.929
# - blandChromatin            1   64.053 78.053
# - maginalAdhesion           1   64.995 78.995
# - bareNuclei                1   75.634 89.634
# - clumpThickness            1   76.942 90.942
# 
# Step:  AIC=68.15
# class ~ clumpThickness + maginalAdhesion + bareNuclei + blandChromatin + 
#   normalNucleoli + mitosis
# 
# Df Deviance    AIC
# <none>                 54.154 68.154
# - mitosis          1   60.177 72.177
# - normalNucleoli   1   60.937 72.937
# - blandChromatin   1   64.056 76.056
# - maginalAdhesion  1   65.022 77.022
# - bareNuclei       1   76.417 88.417
# - clumpThickness   1   77.027 89.027


#========================================================================
# 17.3.1. Classical decision trees
# code listing Creating a classical decision tree with rpart()
library(rpart)
set.seed(1234)
dtree <- rpart(class ~ ., data = df.train, method = "class",
               parms = list(split="information"))
print(dtree)
# n= 489 
# 
# node), split, n, loss, yval, (yprob)
# * denotes terminal node
# 
# 1) root 489 170 benign (0.65235174 0.34764826)  
# 2) sizeUniformity< 2.5 304   8 benign (0.97368421 0.02631579)  
# 4) clumpThickness< 6.5 297   3 benign (0.98989899 0.01010101) *
#   5) clumpThickness>=6.5 7   2 malignant (0.28571429 0.71428571) *
#   3) sizeUniformity>=2.5 185  23 malignant (0.12432432 0.87567568)  
# 6) bareNuclei< 2.5 28  13 benign (0.53571429 0.46428571)  
# 12) sizeUniformity< 3.5 16   1 benign (0.93750000 0.06250000) *
#   13) sizeUniformity>=3.5 12   0 malignant (0.00000000 1.00000000) *
#   7) bareNuclei>=2.5 157   8 malignant (0.05095541 0.94904459) *


dtree$cptable
#           CP nsplit  rel error    xerror       xstd
# 1 0.81764706      0 1.00000000 1.0000000 0.06194645
# 2 0.04117647      1 0.18235294 0.1823529 0.03169642
# 3 0.01764706      3 0.10000000 0.1588235 0.02970979
# 4 0.01000000      4 0.08235294 0.1235294 0.02637116

plotcp(dtree) # figure 17.1

dtree.pruned <- prune(dtree, cp=0.01764706)

library(rpart.plot)
prp(dtree.pruned, type = 2, extra = 104,
    fallen.leaves = TRUE, main="Decision Tree") # figure 17.2

dtree.pred <- predict(dtree.pruned, df.validate, type = "class")
dtree.perf <- table(df.validate$class, dtree.pred,
                    dnn = c("Actual", "Predicted"))
dtree.perf
#            Predicted
# Actual      benign malignant
# benign        129        10
# malignant       4        67

summary(dtree.perf)
# Number of cases in table: 210 
# Number of factors: 2 
# Test for independence of all factors:
#   Chisq = 153.78, df = 1, p-value = 2.585e-35

(129+67)/(129+10+4+67)*100
# [1] 93.33333


# 17.3.2 Conditional inference trees
library(party)
fit.ctree <- ctree(class ~ ., data = df.train)
plot(fit.ctree, main="Conditional Inference Tree")

ctree.pred <- predict(fit.ctree, df.validate, type="response")
ctree.perf <- table(df.validate$class, ctree.pred,
                    dnn=c("Actual", "Predicted"))
ctree.perf
#            Predicted
# Actual      benign malignant
# benign         131         8
# malignant        4        67

summary(ctree.perf)
# Number of cases in table: 210 
# Number of factors: 2 
# Test for independence of all factors:
#   Chisq = 160.72, df = 1, p-value = 7.875e-37

(131+67)/(131+4+67+8)*100
# [1] 94.28571

# 17.4 Random forests
# code listing 17.5. Random forest
library(randomForest)
set.seed(1234)
fit.forest <- randomForest(class ~ ., data = df.train,
                           na.action = na.roughfix,
                           importance=TRUE)
fit.forest
# Call:
# randomForest(formula = class ~ ., data = df.train, importance = TRUE,      na.action = na.roughfix) 
#                Type of random forest: classification
#                      Number of trees: 500
# No. of variables tried at each split: 3
# 
#       OOB estimate of  error rate: 2.66%
# Confusion matrix:
#           benign malignant class.error
# benign       313         6  0.01880878
# malignant      7       163  0.04117647

importance(fit.forest, type = 2)
#                          MeanDecreaseGini
# clumpThickness                  11.382432
# sizeUniformity                  63.037519
# shapeUniformity                 49.027649
# maginalAdhesion                  4.275345
# singleEpithelialCellSize        14.504981
# bareNuclei                      42.736364
# blandChromatin                  22.484755
# normalNucleoli                  11.375285
# mitosis                          1.755382


forest.pred <- predict(fit.forest, df.validate)
forest.perf <- table(df.validate$class, forest.pred,
                     dnn = c("Actual", "Predicted"))
forest.perf
#            Predicted
# Actual      benign malignant
#   benign      128         7
#   malignant     2        68

summary(forest.perf)
# Number of cases in table: 205 
# Number of factors: 2 
# Test for independence of all factors:
#   Chisq = 168.02, df = 1, p-value = 2.004e-38

(128+68)/(128+7+2+68)*100
# [1] 95.60976

# 17.5. Support vector machines
# code listing 17.6. A support vector machine
library(e1071)
set.seed(1234)
fit.svm <- svm(class ~ ., data = df.train)
fit.svm
# Call:
# svm(formula = class ~ ., data = df.train)
# 
# 
# Parameters:
#   SVM-Type:  C-classification 
# SVM-Kernel:  radial 
#       cost:  1 
# 
# Number of Support Vectors:  74


svm.pred <- predict(fit.svm, na.omit(df.validate))
svm.perf <- table(na.omit(df.validate)$class,
                  svm.pred, dnn = c("Actual", "Predicted"))
svm.perf
#            Predicted
# Actual      benign malignant
#   benign       126         9
#   malignant      0        70


# code listing Tuning an RBF support vector machine
set.seed(1234)
tuned <- tune.svm(class ~ ., data = df.train,
                  gamma = 10^(-6:1),
                  cost = 10^(-10:10))
tuned
# Parameter tuning of ‘svm’:
#   
# - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   gamma cost
#    0.01    1
# 
# - best performance: 0.02740302


fit.svm <- svm(class ~ ., data = df.train, gamma = .01, cost = 1)
svm.pred <- predict(fit.svm, na.omit(df.validate))
svm.perf <- table(na.omit(df.validate)$class,
                  svm.pred, dnn = c("Actual", "Predicted"))
svm.perf
#            Predicted
# Actual      benign malignant
#   benign       128         7
#   malignant      0        70

summary(svm.perf)
# Number of cases in table: 205 
# Number of factors: 2 
# Test for independence of all factors:
#   Chisq = 176.7, df = 1, p-value = 2.546e-40

(128+70)/(128+7+70)*100
# [1] 96.58537

# 17.6. Choosing a best predictive solution
# code listing 17.8. Function for assessing binary classification accuracy
performance <- function(table, n=2) {
  if(!all(dim(table) == c(2,2)))
    stop("Must be a 2 x 2 table")
  tn = table[1,1]
  fp = table[1,2]
  fn = table[2,1]
  tp = table[2,2]
  sensitivity = tp/(tp+fn)
  specificity = tn/(tn+fp)
  ppp = tp/(tp+fp)
  npp = tn/(tn+fn)
  hitrate = (tp+tn)/(tp+tn+fp+fn)
  result <- paste("Sensitivity = ", round(sensitivity, n) ,
                  "\nSpecificity = ", round(specificity, n),
                  "\nPositive Predictive Value = ", round(ppp, n),
                  "\nNegative Predictive Value = ", round(npp, n),
                  "\nAccuracy = ", round(hitrate, n), "\n", sep="")
  cat(result)
  }


# code listing 17.9. Performance of breast cancer data classifiers
performance(logit.perf)
# Sensitivity = 0.99
# Specificity = 0.96
# Positive Predictive Value = 0.92
# Negative Predictive Value = 0.99
# Accuracy = 0.97

performance(dtree.perf)
# Sensitivity = 0.94
# Specificity = 0.93
# Positive Predictive Value = 0.87
# Negative Predictive Value = 0.97
# Accuracy = 0.93

performance(ctree.perf)
# Sensitivity = 0.94
# Specificity = 0.94
# Positive Predictive Value = 0.89
# Negative Predictive Value = 0.97
# Accuracy = 0.94

performance(forest.perf)
# Sensitivity = 0.97
# Specificity = 0.95
# Positive Predictive Value = 0.91
# Negative Predictive Value = 0.98
# Accuracy = 0.96

performance(svm.perf)
# Sensitivity = 1
# Specificity = 0.95
# Positive Predictive Value = 0.91
# Negative Predictive Value = 1
# Accuracy = 0.97


#==============================================================
# 17.7 Using the rattle package for data mining
# install.packages("rattle"), detail in the "install_rattle.R"


# The Pima Indians Diabetes dataset is no longer available due to permission restrictions
# loc <- "http://archive.ics.uci.edu/ml/machine-learning-databases/"
# ds <- "pima-indians-diabetes/pima-indians-diabetes.data"

# url <- paste(loc, ds, sep="")
# 
# diabetes <- read.table(url, sep=",", header=FALSE)
# names(diabetes) <- c("npregant", "plasma", "bp", "triceps",
#                      "insulin", "bmi", "pedigree", "age", "class")
# diabetes$class <- factor(diabetes$class, levels=c(0,1),
#                          labels=c("normal", "diabetic"))


library(rattle)
# Loading required package: tibble
# Loading required package: bitops
# Rattle: A free graphical interface for data science with R.
# Version 5.4.0 Copyright (c) 2006-2020 Togaware Pty Ltd.
# Type 'rattle()' to shake, rattle, and roll your data.
# 
# Attaching package: ‘rattle’
# 
# The following object is masked from ‘package:randomForest’:
#   
#   importance


rattle()
# Loading required package: RGtk2
# 2020-11-14 11:03:35.908 rsession[43141:820420] *** WARNING: Method userSpaceScaleFactor 
# in class NSView is deprecated on 10.7 and later. It should not be used in new 
# applications. Use convertRectToBacking: instead.
#
# Attaching package: ‘zoo’

# The following objects are masked from ‘package:base’:
  
#   as.Date, as.Date.numeric

# export a Decision Tree model with pmml package.
install.packages('pmml') 
```
