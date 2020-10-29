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

##### 8.2.2. Simple linear regression

- Figure 8.1. Scatter plot with regression line for weight predicted from height

- with `lm()` and `abline()`
  
  ![listing1](listing1.png)

##### 8.2.3. Polynomial regression

- Figure 8.2. Quadratic regression for weight predicted by height
  
  ![listing2](listing2.png)

- Figure 8.3. Scatter plot of height by weight, with linear and smoothed fits, and marginal box plots

- with **scatterplot()** in the "car" package

  ![listingn](listingn.png)

##### 8.2.4. Multiple linear regression

- Figure 8.4. Scatter plot matrix of dependent and independent variables for the states data, including linear and smoothed fits, and marginal distributions (kernel density plots and rug plots)
  
  ![listing3](listing3.png)

##### 8.2.5. Multiple linear regression with interactions

- Figure 8.5. Interaction plot for hp*wt. This plot displays the relationship between mpg and hp at 3 values of wt.

  ![effplot](effplot.png)

#### 8.3. Regression diagnostics

##### 8.3.1. A typical approach

- Figure 8.6. Diagnostic plots for the regression of weight on height

  ![typapp](typ821.png)

- Figure 8.7. Diagnostic plots for the regression of weight on height and height-squared

  ![typ8311](typ8311.png)

- Figure 8.8. Diagnostic plots for the regression of murder rate on state characteristics

  ![typ8312](typ8312.png)

##### 8.3.2. An enhanced approach

- library(car)
  
  ![tabenc](tabenc.png)

- Normality
  
  - Figure 8.9. Q-Q plot for studentized residuals
  
  ![fig89](fig89.png)

  - Figure 8.10. Distribution of studentized residuals using the residplot() function
  
  - with residplot()

  ![fig810](fig810.png)

- Independence of Errors
  
  - with durbinWatsonTest()

- Linearity

  - Figure 8.11. Component plus residual plots for the regression of murder rate on state characteristics

  ![fig811](fig811.png)

- Homoscedasticity

  - Figure 8.12. Spread-level plot for assessing constant error variance

  ![fig812](fig812.png)

##### 8.3.3. Global validation of linear model assumption

- gvlma() in the gvlma package.
  
##### 8.3.4. Multicollinearity

- vif()

#### 8.4. Unusual observations

##### 8.4.1. Outliers

- outlierTest()
  
##### 8.4.2. High leverage points

- Figure 8.13. Index plot of hat values for assessing observations with high leverage

  ![fig813](fig813.png)

##### 8.4.3. Influential observations

- Figure 8.14. Cook’s D plot for identifying influential observations

  ![fig814](fig814.png)

- Figure 8.15. Added-variable plots for assessing the impact of influential observations

  ![fig815](fig815.png)

- Figure 8.16. Influence plot.

  - States above `+2` or below `–2` on the vertical axis are considered **outliers**.
  
  - States above `0.2` or `0.3` on the horizontal axis **have high leverage** (unusual combinations of predictor values).
  
  - **Circle size** is proportional to influence. Observations depicted by large circles may have disproportionate influence on the parameters estimates of the model.

  ![fig816](fig816.png)

#### 8.5. Corrective measures

##### 8.5.1. Deleting observations

##### 8.5.2. Transforming variables

- with powerTransform() and boxTidwell()

##### 8.5.3. Adding or deleting variables

- sqrt(vif) > 2

##### 8.5.4. Trying a different approach

#### 8.6. Selecting the “best” regression model

- The selection of a final regression model always involves a compromise between predictive accuracy (a model that fits the data as well as possible) and parsimony (a simple and replicable model).

##### 8.6.1. Comparing models

- with anova() and AIC()

##### 8.6.2. Variable selection

- selecting a final set of predictor variables from a larger pool of candidate variables

  - Stepwise Regression
    - stepAIC() in the MASS package
  
  - all-subsets regression
    - regsubsets() in the leaps package
    - Figure 8.17. Best four models for each subset size based on Adjusted R-square
    ![fig817](fig817.png)

    - Figure 8.18. Best four models for each subset size based on the Mallows Cp statistic
    ![fig818](fig818.png)

#### 8.7. Taking the analysis further

##### 8.7.1. Cross-validation

- crossval() in the bootstrap package

##### 8.7.2. Relative importance

- Figure 8.19. Bar plot of relative weights for the states multiple regression problem
  
  ![fig819](fig819.png)
