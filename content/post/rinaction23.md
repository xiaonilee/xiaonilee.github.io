---
title: "Chapter 23. Advanced graphics with the lattice package"
date: 2020-11-23
lastmod: 2020-11-23
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

Notebook of Reading Books: R in Action_Chapter 23.

<!--more-->

## This chapter covers

- An introduction to the lattice package

- Grouping and conditioning

- Adding information with panel functions

- Customizing a lattice graph’s appearance

### 23.1. The lattice package

- Figure 23.1. Trellis graph of singer heights by voice part.

  ![fig231](fig231.png)

- Table 23.1. Graph types and corresponding functions in the ***lattice*** package.

  ![tab231](tab231.png)

### 23.2. Conditioning variables

- Typically, conditioning variables are factors. 

- If It's a continuous variable:
  - `cut()`
  - **myshingle**

```r
myshingle <- equal.count(x, number=n, overlap=proportion)
```

- Figure 23.2. Trellis plot of miles per gallon vs. car weight conditioned on engine displacement.
  - engine displacement is a continuous variable
    - it has been converted to three non-overlapping shingles with equal numbers of observations.

  ![fig232](fig232.png)

### 23.3. Panel functions

- Figure 23.3. Trellis plot of miles per gallon vs. car weight conditioned on engine displacement. 
  - A custom panel function has been used to add regression lines, rug plots, and grid lines.

  ![fig233](fig233.png)

- Figure 23.4. Trellis graph of miles per gallon vs. engine displacement conditioned on transmission type. 
  - Smoothed lines (loess), grids, and group mean levels have been added.

  ![fig234](fig234.png)

### 23.4. Grouping variables

Attach is the [Script](chapter23.R) of chapter23.

Show me the code <i class="far fa-hand-pointer"></i>

```r

```
