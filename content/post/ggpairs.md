---
title: "Review raw data with ggpairs"
date: 2021-03-21
lastmod: 2021-03-21
draft: false
tags: ["Bioinformatics", "R", "R package"]
categories: ["Bioinformatics", "R", "R package"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

It's very necessary to have a look at the data source before trying to analyze it.

In this blog, I will share a powerful function `ggpairs` to implement the purpose. 


<!--more-->

### Package

- **GGaly**

- **ggplot2**

- **reshape**

### Data Source

- `tips` from `reshape` package


### Load data

```r
library(GGally)
library(ggplot2)
data(tips, package = "reshape")
head(tips)
colnames(tips)
dim(tips)
```

![1](1.png)

### For each column

```r
ggpairs(tips, mapping = aes(color = sex))
```

![2](2.png)

### Select interested columns

```r
ggpairs(tips, mapping = aes(color = sex), 
        columns = c("total_bill", "time", "tip"), 
        columnLabels = c("Total_Bill(Continuous_variable)", "Time(discrete_variable)", "Tip(Continuous_variable)"),
        title = "Select_column")
```

![3](3.png)
