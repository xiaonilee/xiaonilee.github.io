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


<!--more-->

## Chapter1. Introduction To R

- This chapter covers
  
  - Installing R
  
  - Understanding the R language
  
  - Running programs

- In the present chapter, a tremendously `versatile` and `impressive` code need to be marked, to remove most objects from the working environment.

```r
rm(list = ls())
options(stringsAsFactors = F)
```

- To summary, R can
  
  - Access data from a wide range of sources
  
  - Merge the pieces of data together
  
  - Clean and annotate them
  
  - Analyze them with the latest methods
  
  - Present the findings in meaningful and graphically appealing ways
  
  - Incorporate the results into attractive reports that can be distributed to stakeholders and the public.

Attach is the [Script](chapter1.R) of chapter1.
