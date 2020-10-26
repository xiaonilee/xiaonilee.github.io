---
title: "Reading Books : R in Action"
date: 2020-10-19
lastmod: 2020-10-25
draft: false
tags: ["R", "Bioinformatics", "Book"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Notebook of Reading Books: R in Action.

<!--more-->

## Chapter1. Introduction To R

### Chapter 1 covers
  
- Installing R
  
- Understanding the R language
  
- Running programs

- In the present chapter, a tremendously `versatile` and `impressive` code need to be marked, which is used to remove most objects from the working environment.

```r
rm(list = ls())
options(stringsAsFactors = F)
```

### To summary, R can
  
- Access data from a wide range of sources;
  
- Merge the pieces of data together;
  
- Clean and annotate them;
  
- Analyze them with the latest methods;
  
- Present the findings in meaningful and graphically appealing ways;
  
- Incorporate the results into attractive reports that can be distributed to stakeholders and the public.

Attach is the [Script](chapter1.R) of chapter1.

## Chapter2. Creating Dataset

### Chapter 2 covers
  
- Exploring R data structures
  
- Using data entry
  
- Importing data
  
- Annotating datasets

- In summary, this chapter describes
  
  - The various structures that R provides for holding data;
  
  - The many methods available for importing data from both keyboard and external sources.

Attach is the [Script](chapter2.R) of chapter2.

## Chapter3. Getting Started with Graphs

### Chapter 3 covers
  
- Creating and saving graphs
  
- Customizing symbols, lines, colors, and axes
  
- Annotating with text and titles#   Controlling a graph’s dimensions
  
- Combining multiple graphs into one

**Figure for Comparing Drug A and Drug B response by dose with code Listing 3.3.**

![chapter3_3_3](chapter3_3_3.png)
  
Attach is the [Script](chapter3.R) of chapter3.

## Chapter4. Basic Data Management

### Chapter 4 covers
  
- Manipulating dates and missing values
  
- Understanding data type conversions
  
- Creating and recoding variables
  
- Sorting, merging, and subsetting datasets
  
- Selecting and dropping variables

Attach is the [Script](chapter4.R) of chapter4.

## Chapter5. Advanced data management

### Chapter 5 covers
  
- Mathematical and statistical functions
  
- Character functions
  
- Looping and conditional execution
  
- User-written functions
  
- Ways to aggregate and reshape data

Attach is the [Script](chapter5.R) of chapter5.

## Chapter6. Basic graphs

### Chapter 6 covers
  
- Bar, box, and dot plots
  
- Pie and fan charts
  
- Histograms and kernel density plots

**Figure for Pie charts with code Listing 6.5.**

![piecharts](piecharts.png)

**Figure for histogram charts with code Listing 6.6.**

- The surrounding box is produced by the `box()` function.

![hist6.6](hist6.6.png)

**Figure for density plots with code Listing 6.7 and 6.8.**

- The `locator(1)` option indicates that you’ll place the legend interactively by clicking on the graph where you want the legend to appear.

![density](densityplot.png)

**Figure for box plots with code listing 6.9.**

![boxplot](boxplot1.png)

**Figure for violin plots with code listing 6.10.**

- A `violin plot` is a **combination** of a ***box plot*** and a ***kernel density plot***.

![violinplot](violinplot.png)

**Figure for dot charts with code listing 6.11.**

- Dot plot of mpg for car models **grouped**, **sorted**, and **colored**.

![dotchart](dotchart.png)
