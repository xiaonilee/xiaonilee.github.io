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

Attach is the [Script](chapter18.R) of chapter18.

Show me the code <i class="far fa-hand-pointer"></i>

```r

```
