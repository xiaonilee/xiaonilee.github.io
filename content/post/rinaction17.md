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

- rattle()

Attach is the [Script](chapter17.R) of chapter17.

Show me the code <i class="far fa-hand-pointer"></i>

```r

```
