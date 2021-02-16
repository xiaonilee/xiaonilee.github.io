---
title: "The way to use R in Google Colab"
date: 2021-02-16
lastmod: 2021-02-16
draft: false
tags: ["R package", "Bioinformatics"]
categories: ["R package", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

**Colab** is an interactive notebook provided by Google (primarily) for writing and running `Python` through a *browser*. 

Then, how to use R in Google Colab?

<!--more-->

### There are two ways to run R in Colab

- The first way is to use the `rpy2` package in the Python runtime. This method allows you to execute R and Python syntax together.

- The second way is to actually start the notebook in the R runtime.

### To use R and Python together in Colab

- Open google browser.

- Create a new [notebook](https://colab.research.google.com/#create=true).

- Run rmagic by executing this command `%load_ext rpy2.ipython`.

- After that, every time you want to use R, add `%%R` in the beginning of each cell.

  - Start rmagic by executing this in a cell:

  ```python
  %load_ext rpy2.ipython
  ```

  - Use %%R to execute cell magic

  ```python
  %%R
  x <- seq(0, 2*pi, length.out=50)
  x
  ```

  - copy R variable to Python:

  ```python
  x = %R x
  ```

### Start to use R in Colab

- Go to this [URL](https://colab.research.google.com/#create=true&language=r), or this short [URL](https://colab.to/r).

- Confirm that you are in the R runtime by going to the “`Runtime`” settings, and select “`Change runtime type`”.

- Check and print out the R version

```r
R.version.string
```

- Check packages available

```r
print(installed.packages())
```

### As an example, to install **Package esquisse**

- `Package esquisse`
  
  - Explore and Visualize Your Data Interactively with the ggplot2 package.
  
  - A 'shiny' gadget. 
  
  - To create 'ggplot2' charts interactively with drag-and-drop to map your variables. 
  
  - You can quickly visualize your data accordingly to their type, export to 'PNG' or 'PowerPoint', and retrieve the code to reproduce the chart.

    ```r
    install.packages("esquisse")
    library(esquisse)
    ```
