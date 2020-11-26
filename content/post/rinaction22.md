---
title: "Chapter 22. Creating dynamic reports"
date: 2020-11-22
lastmod: 2020-11-23
draft: false
tags: ["R", "R in Action", "Bioinformatics", "Book"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

Notebook of Reading Books: R in Action_Chapter 22.

<!--more-->

## This chapter covers

- Publishing results to the web

- Incorporating R results into Microsoft Word or Open Document reports

- Creating dynamic reports, where changing the data changes the report

- Creating publication quality documents with R, Markdown, and LaTeX

### 22.1. A template approach to reports 

- Four types of template
  - R Markdown
  - ODT
  - DOCX
  - LaTex

### 22.2. Creating dynamic reports with R and Markdown

- Step by Step
  - install.packages("rmarkdown) for older versions of RStudio
  - install.packages("xtable"), which can attractively formats data frames and matrices for inclusion in reports. 
  - install MacTeX software to generate PDF documentation

- To generate a [html](women.html) document for the [women.Rmd](women.Rmd) template.

  ![womenhtml](womenhtml.png)

- To generate a nicely formatted [PDF](womenPDF.pdf) document for the [womenPDF.Rmd](womenPDF.Rmd) template.

  ![womenpdf.pdf](womenpdf.png)

- To generate an attractive [word](womenWord.docx) document for the [womenWord.Rmd](womenWord.Rmd) template.

  ![womenword](womenword.png)

### 22.3. Creating dynamic reports with R and LaTeX

- LaTeX is powerful, But it requires significant study to use effectively and creates documents in formats (PDF, DVI, PS) that can’t be edited.

- With the text file [**drugs.Rnw**](drugs.Rnw):
  - Processed through the `knit()` function.
    - Resulting in the output file: [drugs.tex](drugs.tex) and [figure](figure/unnamed-chunk-4-1.pdf).

      ![fig](fig.png)

  - Processed through the `knit2pdf()` function.
    - Resulting in a typeset PDF document ([drugs.pdf](drugs.pdf)).

      ![drugpdf](drugpdf.png)

### 22.4. Creating dynamic reports with R and Open Document

- In the `odfWeave()` package.

### 22.5. Creating dynamic reports with R and Microsoft Word

- It's only work on Windows platforms.
- In the `R2wd()` package.

Attach is the [Script](chapter22.R) of chapter22.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)


# 22.2. Creating dynamic reports with R and Markdown
# code listing 22.1. The document of women.Rmd: a Markdown template with embedded R code

# render the file
# setwd() to current path
library(rmarkdown)
render("women.Rmd", "html_document")

# render the file
library(rmarkdown)
# A new version of TeX Live has been released, so
# tinytex::reinstall_tinytex()
render("womenPDF.Rmd", "pdf_document")

# render the file
library(rmarkdown)
render("womenWord.Rmd", "word_document")


# 22.3. Creating dynamic reports with R and LaTeX
# code listing 22.2. drugs.Rnw: a sample LaTeX template with embedded R code

library(knitr)
knit("drugs.Rnw")

library(knitr)
knit2pdf("drugs.Rnw")
```
