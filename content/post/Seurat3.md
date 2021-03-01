---
title: "Analysis, visualization, and integration of spatial datasets with Seurat"
date: 2021-02-28
lastmod: 2021-03-01
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

how to use Seurat to analyze spatially-resolved RNA-seq data?

Herein, the tutorial will cover these tasks:

- Normalization

- Dimensional reduction and clustering

- Detecting spatially-variable features

- Interactive visualization

- Integration with single-cell RNA-seq data

- Working with multiple slices

<!--more-->

### Material and Methods

- **Dataset**: sagital mouse brain slices generated using the Visium v1 chemistry

### Load dataset

```r
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#BiocManager::install("rhdf5")
#BiocManager::install("rhdf5r")
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
```

### Data preprocessing
- Raw data performance
  
    ```r
    plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
    wrap_plots(plot1, plot2)
    ggsave("./fig1.png")
    ```

    ![fig1](fig1.png)

- The variance in molecular counts / spot can be substantial for spatial datasets, particularly if there are differences in cell density across the tissue. 

- Normalize the data in order to account for variance in sequencing depth across data points

    ```r
    # Normalizes the data, detects high-variance features, and stores the data in the SCT assay.
    brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
    ```

### Gene expression visualization

- `SpatialFeaturePlot()` function can overlay molecular data on top of tissue histology. 

    ```r
    SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
    ggsave("./fig2.png")
    ```

    ![fig2]("./fig2.png")

- Adjust the size of the spots (and their transparency) to improve the visualization of the histology image.
  
    ```r
    p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)

    # to downweight the transparency of points with lower expression with alpha
    p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
    p1 + p2
    ggsave("./fig3.png")
    ```

    ![fig3](fig3.png)

### Dimensionality reduction, clustering, and visualization

- With the same workflow to run PCA and clustering

```r
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
```

- Visualize the results

    ```r
    #in UMAP space  with DimPlot()
    p1 <- DimPlot(brain, reduction = "umap", label = TRUE)

    #overlaid on the image with SpatialDimPlot()
    p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
    p1 + p2
    ```

    ![fig4](fig4.png)

- Visualize special voxel belongs to which cluster

    ```r
    SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 
        5, 8)), facet.highlight = TRUE, ncol = 3)
    ggsave("./fig5.png")
    ```

    ![fig5](fig5.png)

### Interactive plotting

- Both `SpatialDimPlot()` and `SpatialFeaturePlot()` now have an interactive parameter, that when set to `TRUE`, open up the Rstudio viewer pane with an interactive Shiny plot.

    ```r
    #with SpatialDimPlot()
    SpatialDimPlot(brain, interactive = TRUE)

    #with SpatialFeaturePlot()
    SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

    #with LinkedDimPlot()
    LinkedDimPlot(brain)
    ggsave("./fig6.png")
    ```

    ![fig6](fig6.png)

### Identification of Spatially Variable Features

- Seurat offers two workflows to identify molecular features that correlate with spatial location within a tissue
  - Based on pre-annotated anatomical regions within the tissue, which may be determined either from unsupervised clustering or prior knowledge.
  
    ```r
    de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
    SpatialFeaturePlot(object = brain, 
                    features = rownames(de_markers)[1:3], 
                    alpha = c(0.1, 1), ncol = 3)
    ggsave("./fig7.png")
    ```

    ![fig7](fig7.png)

  - An alternative approach, implemented in `FindSpatiallyVariables()`.
  
    ```r
    brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", 
                                        features = VariableFeatures(brain)[1:1000], 
                                        selection.method = "markvariogram")
    ggsave("./fig8.png")
    ```

- Visualize the expression of the top 6 features identified by this measure.

    ```r
    top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
    SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
    ggsave("./fig9.png")
    ```

### Subset out anatomical regions

```r
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 |
# image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
ggsave("./fig10.png")
```

![fig10](fig10.png)

### Integration with single-cell data

- Download [data](https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1)

    ```r
    allen_reference <- readRDS("/Users/xiaonili/Downloads/allen_cortex.rds")
    # note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells
    # this speeds up SCTransform dramatically with no loss in performance
    library(dplyr)
    allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)

    # After subsetting, we renormalize cortex
    cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
    # the annotation is stored in the 'subclass' column of object metadata
    DimPlot(allen_reference, group.by = "subclass", label = TRUE)
    ggsave("./fig11.png")
    ```

    ![fig11](fig11.png)

- Get prediction scores for each spot for each class.

```r
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE)
ggsave("./fig12.png")
```

![fig12](fig12.png)

- Based on these prediction scores, to predict `cell types` whose location is spatially restricted.

    ```r
    cortex <- FindSpatiallyVariableFeatures(cortex, 
                                            assay = "predictions", 
                                            selection.method = "markvariogram", 
                                            features = rownames(cortex), 
                                            r.metric = 5, slot = "data")
    top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
    SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
    ggsave("./fig13.png")
    ```

    ![fig13](fig13.png)

- Finally, we show that our integrative procedure is capable of recovering the known spatial localization patterns of both neuronal and non-neuronal subsets, including laminar excitatory, layer-1 astrocytes, and the cortical grey matter.

    ```r
    SpatialFeaturePlot(cortex, 
                    features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", 
                                            "L6b", "Oligo"), 
                    pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
    ggsave("./fig14.png")
    ```

    ![fig14](fig14.png)


### Working with multiple slices in Seurat

- Load datasets

```r
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
```

- Merge and work multiple slices

    ```r
    brain.merge <- merge(brain, brain2)
    ```

    ```r
    DefaultAssay(brain.merge) <- "SCT"
    VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
    brain.merge <- RunPCA(brain.merge, verbose = FALSE)
    brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
    brain.merge <- FindClusters(brain.merge, verbose = FALSE)
    brain.merge <- RunUMAP(brain.merge, dims = 1:30)
    ```

    ```r
    DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
    ggsave("./fig15.png")
    ```

    ![fig15](fig15.png)

    ```r
    SpatialDimPlot(brain.merge)
    ggsave("./fig16.png")
    ```

    ![fig16](fig16.png)

    ```r
    SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))
    ggsave("./fig17.png")
    ```

    ![fig17](fig17.png)

### Slide-seq

- Load datasets

    ```r
    InstallData("ssHippo")
    slide.seq <- LoadData("ssHippo")
    ```

- Data preprocessing

    ```r
    plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
    slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
    plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
    wrap_plots(plot1, plot2)
    ggsave(""./fig18.png)
    ```

    ![fig18](fig18.png)

- Normalize

    ```r
    slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
    slide.seq <- RunPCA(slide.seq)
    slide.seq <- RunUMAP(slide.seq, dims = 1:30)
    slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
    slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
    ```

- Visualize

    ```r
    plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
    plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
    plot1 + plot2
    ggsave("./fig19.png")
    ```


- Alternately,

    ```r
    SpatialDimPlot(slide.seq, 
                cells.highlight = CellsByIdentities(object = slide.seq, 
                                                    idents = c(1, 6, 13)), 
                facet.highlight = TRUE)
    ggsave("./fig20.png")     
    ```





