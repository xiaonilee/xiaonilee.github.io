---
title: "Using Seurat with multimodal data"
date: 2021-02-26
lastmod: 2021-02-26
draft: false
tags: ["Bioinformatics", "R package", "RNA Sequence", "PCA", "KNN"]
categories: ["Bioinformatics", "R package", "RNA Sequence", "PCA", "KNN"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

The ability to make simultaneous measurements of multiple data types from the same cell, known as multimodal analysis, represents a new and exciting frontier for single-cell genomics.  
Seurat4 to enable for the seamless storage, analysis, and exploration of diverse multimodal single-cell datasets.


Herein, I will follow official Tutorial for using Seurat with multimodal data step by step.

<!--more-->

### Metarial and Methods

- **Dataset**: 8,617 cord blood mononuclear cells (CBMCs)


### Load in the data

- Link for the [RNA UMI matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz) download.

- Link for the [ADT UMI matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz) download.

```r
library(Seurat)
library(ggplot2)
library(patchwork)
```

```r
# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls
# for the protein measurements. For this reason, the gene expression matrix has HUMAN_ or MOUSE_
# appended to the beginning of each gene.
cbmc.rna <- as.sparse(read.csv(file = "./GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
    header = TRUE, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "./GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", 
    header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))
```

![fig1](1.png)
