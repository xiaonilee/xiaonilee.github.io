---
title: "Chapter 1. Introduction To R"
date: 2020-10-19
lastmod: 2020-10-28
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

Notebook of Reading Books: R in Action_Chapter 1.

<!--more-->

### This chapter covers
  
- Installing R
  
- Understanding the R language
  
- Running programs

### To summary, R can
  
- **Access** data from a wide range of sources;
  
- **Merge** the pieces of data together;
  
- **Clean** and **annotate** data;
  
- **Analyze** cleaned data with the latest methods;
  
- Present the findings in **meaningful** and **graphically appealing** ways;
  
- Incorporate the results into attractive reports that can be distributed to stakeholders and the public.

Attach is the [Script](chapter1.R) [of](script1.R) [chapter](script2.R)[1](script3.R).

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
# a tremendously versatile and impressive code
options(stringsAsFactors = F)

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
