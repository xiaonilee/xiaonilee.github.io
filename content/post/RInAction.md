---
title: "Reading Books : R in Action"
date: 2020-10-19
lastmod: 2020-10-27
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

Notebook of Reading Books: R in Action.

<!--more-->

## Chapter1. Introduction To R

### Chapter 1 covers
  
- Installing R
  
- Understanding the R language
  
- Running programs

- In the present chapter, a tremendously **versatile** and **impressive** code need to be marked, which is used to remove most objects from the working environment.

```r
rm(list = ls())
options(stringsAsFactors = F)
```

### To summary, R can
  
- **Access** data from a wide range of sources;
  
- **Merge** the pieces of data together;
  
- **Clean** and **annotate** data;
  
- **Analyze** cleaned data with the latest methods;
  
- Present the findings in **meaningful** and **graphically appealing** ways;
  
- Incorporate the results into attractive reports that can be distributed to stakeholders and the public.

Attach is the [Script](chapter1.R) of chapter1.

Show me the code

```r
# Section 1.3.2 Getting help
help.start()

## Section 1.3.3 Functions for managing the R workspace
getwd()
setwd("mydirectory")
ls()            # List the objects in the current workspace.
rm(objectlist)  # Remove (delete) one or more objects.
options()       # View or set current options.
history()       # Display your last # commands (default = 25).
save()
load()
q()

## Section 1.3.4 Input and Output

# Input
source("filename")

source("script1.R") # Working through an example.

# Output
sink("myoutput", append=TRUE, split=TRUE)
pdf("mygraphs.pdf")
source("script2.R")

sink()
dev.off()
source("script3.R")
```

## Chapter2. Creating Dataset

### Chapter 2 covers
  
- Exploring R data structures
  
- Using data entry
  
- Importing data
  
- Annotating datasets

- In summary, this chapter describes
  
  - The various structures that R provides for holding data;

  ![datastructures](datastructures.png)
  
  - The many methods available for importing data from both keyboard and external sources.

Attach is the [Script](chapter2.R) of chapter2.

## Chapter3. Getting Started with Graphs

### Chapter 3 covers
  
- Creating and saving graphs
  
- Customizing symbols, lines, colors, and axes
  
- Annotating with text and titles
  
  - `text()` and `mtext()`

- Controlling a graph’s dimensions

  - `par()`
  
- Combining multiple graphs into one

  - `help("layout")`

### A Example about Graphical Parameters

- **Figure for Comparing Drug A and Drug B response by dose with code Listing 3.3.**

- ***par()***

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

### Notes about some useful functions

- 5.2.2. Statistical functions

  - scale()

  - tansform()

- 5.2.3. Probability functions

  - runif()

  - set.sead()

  - mvrnorm()

- 5.2.4. Character functions

  - nchar()

  - substr()

  - grep()

  - sub()

  - strsplit()

  - paste()

  - toupper() and tolower()

  - seq(from, to, by)

  - rep()

  - cut()

  - pretty()

  - cat()

- 5.2.6. Applying functions to matrices and data frames

  - apply(x, margin, FUN, ...)
  
- 5.4. Control flow

  - for (var in seq) statement

  - if (cond) statement [else statement2]

  - ifelse(cond, statement1, statement2)

  - switch(expr, ...)

- 5.6.2. Aggregating data

  - aggregate(x, by, FUN)

- 5.6.3. The reshape package, versatility

  - melt()

  - cast() and dcast()

Attach is the [Script](chapter5.R) of chapter5.

## Chapter6. Basic graphs

### Chapter 6 covers
  
- Bar, box, and dot plots
  
- Pie and fan charts
  
- Histograms and kernel density plots

### Note about pie charts

- **Figure for Pie charts with code Listing 6.5.**

![piecharts](piecharts.png)

### Note about histogram charts

- **Figure for histogram charts with code Listing 6.6.**

- The surrounding box is produced by the ***box()*** function.

![hist6.6](hist6.6.png)

### Note about density plots

- **Figure for density plots with code Listing 6.7 and 6.8.**

- The ***locator(1)*** option indicates that you’ll place the legend interactively by clicking on the graph where you want the legend to appear.

![density](densityplot.png)

### Note about box plots

- **Figure for box plots with code listing 6.9.**

![boxplot](boxplot1.png)

### Note about violin plots

- **Figure for violin plots with code listing 6.10.**

- A `violin plot` is a **combination** of a ***box plot*** and a ***kernel density plot***.

![violinplot](violinplot.png)

### Note about dot charts

- **Figure for dot charts with code listing 6.11.**

- Dot plot of mpg for car models **grouped**, **sorted**, and **colored**.

![dotchart](dotchart.png)

Attach is the [Script](chapter6.R) of chapter6.

## Chapter7. Basic statistics

### Chapter 7 covers

- Descriptive statistics
  
- Frequency and contingency tables
  
- Correlations and covariances
  
- t-tests
  
- Nonparametric statistics

#### In summary

- 7.1.1 Descriptive statistics

  - summary()

  - sapply()

  - describe() in the Hmisc package()

  - stat.desc() in the pastecs package

  - describe() in the psych package
  
- 7.1.2 Descriptive statistics by group
  
  - aggregate()
  
  - by()
  
  - summaryBy() in the doBy package

  - describe.by() in the psych package

- 7.2.1 Generating frequency tables

  - table()

  - prop.table()

  - xtabs()

  - CrossTable()

  - ftable()

- 7.2.2 Tests of independence
  
  - chisq.test()

  - fisher.test()

  - mantelhaen.test()

- 7.2.3 Measures of association

  - assocstats()

- 7.3 Correlations

  - cov()

  - cor()

  - pcor() in the ggm package

  - cor.test()

  - corr.test in the psych package

- 7.4 t-tests

  - t.test()

- 7.5 Nonparametric tests of group differences

  - wilcox.test()

  - kruskal.test()
  
Attach is the [Script](chapter7.R) of chapter7.
