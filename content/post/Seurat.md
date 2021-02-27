---
title: "Seurat - Guided Clustering Tutorial"
date: 2021-02-25
lastmod: 2021-02-27
draft: false
tags: ["Bioinformatics", "R", "RNA Sequence", "PCA", "KNN"]
categories: ["Bioinformatics", "R", "RNA Sequence", "PCA", "KNN"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

In this article, I will follow the official Tutorial to do `clustering` using Seurat step by step.


<!--more-->


![The finnal goal is the figure](dimplot4.png)


### Metarial and Methods

- **Dataset**: a dataset of 2700 Peripheral Blood Mononuclear Cells freely available from 10X Genomics

- **Flatform**: Illumina NextSeq 500
### Setup the Seurat Object

```r
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
```

```r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size
```

  ![fig1](1.png)

### Standard pre-processing workflow

- QC and selecting cells for further analysis

```r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("./vlnplot.png")
```

![vlnplot](vlnplot.png)

```r
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("./plot12.png")
```

![plot12](plot12.png)

```r
#filter cells that have unique feature counts over 2,500 or less than 200 & have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

### Normalizing the data

```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
```

### Identification of highly variable features (feature selection)

```r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

![plot1212](plot1212.png)

### Scaling the data

```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

### Perform linear dimensional reduction

```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("./dimReduction.png")
```

![dimreduction](dimReduction.png)

```r
DimPlot(pbmc, reduction = "pca")
ggsave("./dimplot.png")
```

![dimplot](dimplot.png)

```r
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ggsave("./dimheatmap1.png")
```

![dimheatmap1](dimheatmap1.png)

```r
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("./dimheatmap2.png")
```

![dimheatmap2](dimheatmap2.png)


### Determine the ‘dimensionality’ of the dataset

```r
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#Comparing the distribution of p-values for each PC with a uniform distribution. 
JackStrawPlot(pbmc, dims = 1:15)
ggsave("./jackstrawplot.png")
```

![jackstrawplot](jackstrawplot.png)

```r
# An alternative heuristic method generates an ‘Elbow plot’
ElbowPlot(pbmc)
ggsave("./elbowplot.png")
```

![elbowplot](elbowplot.png)

### Cluster the cells

```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```

![2](2.png)

### Run non-linear dimensional reduction (UMAP/tSNE)

```r
# reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
ggsave("./dimplot2.png")

#saveRDS(pbmc, file = "./output/pbmc_tutorial.rds")
```

![dimplot2](dimplot2.png)

### Finding differentially expressed features (cluster biomarkers)

```r
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

![3](3.png)

```r
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
ggsave("./vlnplot2.png")
```

![vlnplot2](vlnplot2.png)

```r
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
ggsave("./vlnplot3.png")
```

![vlnplot3](vlnplot3.png)

```r
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
ggsave("./featureplot.png")
```

![featureplot](featureplot.png)

```r
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("./doheatmap.png")
```

![doheatmap](doheatmap.png)


### Assigning cell type identity to clusters

```r
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("./dimplot4.png")
#saveRDS(pbmc, file = "./output/pbmc3k_final.rds")
```

![dimplot4](dimplot4.png)
