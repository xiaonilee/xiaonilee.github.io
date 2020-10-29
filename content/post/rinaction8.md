---
title: "Chapter 8. Regression"
date: 2020-10-29
lastmod: 2020-10-29
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

Notebook of Reading Books: R in Action_Chapter 8.

<!--more-->

### This chapter covers

- Fitting and interpreting linear models

- Evaluating model assumptions

- Selecting among competing models
  
#### 8.1. The many faces of regression

#### 8.2. OLS regression

In OLS regression, a quantitative dependent variable is predicted from a weighted sum of predictor variables, where the weights are parameters estimated from the data.

$$ \begin{aligned} \hat{Y}_i & = \beta_0 + \beta_1 X_1i\ + ... + \beta_kX_ki\ & i = 1...n \end{aligned} $$

##### 8.2.1. Fitting regression models with lm()

![tabfun](tabfun.png)

#### 8.2.2. Simple linear regression

- with `lm()` and `abline()`
  
  ![listing1](listing1.png)

#### 8.2.3. Polynomial regression

  ![listing2](listing2.png)

- with **scatterplot()** in the "car" package

  ![listingn](listingn.png)

#### 8.2.4. Multiple linear regression

  ![listing3](listing3.png)

#### 8.2.5. Multiple linear regression with interactions

- Interaction plot for hp*wt. This plot displays the relationship between mpg and hp at 3 values of wt.

  ![effplot](effplot.png)

#### 8.3. Regression diagnostics
