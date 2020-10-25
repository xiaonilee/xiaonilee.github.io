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

- This chapter covers
  
  - Installing R
  
  - Understanding the R language
  
  - Running programs

- In the present chapter, a tremendously `versatile` and `impressive` code need to be marked, which is used to remove most objects from the working environment.

```r
rm(list = ls())
options(stringsAsFactors = F)
```

- To summary, R can
  
  - Access data from a wide range of sources;
  
  - Merge the pieces of data together;
  
  - Clean and annotate them;
  
  - Analyze them with the latest methods;
  
  - Present the findings in meaningful and graphically appealing ways;
  
  - Incorporate the results into attractive reports that can be distributed to stakeholders and the public.

Attach is the [Script](chapter1.R) of chapter1.

## Chapter2. Creating Dataset

- This chapter covers
  
  - Exploring R data structures
  
  - Using data entry
  
  - Importing data
  
  - Annotating datasets

- In summary, this chapter describes
  
  - The various structures that R provides for holding data;
  
  - The many methods available for importing data from both keyboard and external sources.

Attach is the [Script](chapter2.R) of chapter2.

## Chapter3. Getting Started with Graphs

- This chapter covers
  
  - Creating and saving graphs
  
  - Customizing symbols, lines, colors, and axes
  
  - Annotating with text and titles#   Controlling a graph’s dimensions
  
  - Combining multiple graphs into one
  
Attach is the [Script](chapter3.R) of chapter3.

