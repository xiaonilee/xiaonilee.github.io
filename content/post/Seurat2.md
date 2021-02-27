---
title: "Using Seurat with multimodal data"
date: 2021-02-26
lastmod: 2021-02-27
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

The ability to make simultaneous measurements of multiple data types from the same cell, known as multimodal analysis, represents a new and exciting frontier for ***single-cell genomics***.  `Seurat4` to enable for the seamless storage, analysis, and exploration of diverse multimodal single-cell datasets.


Herein, I will follow the official Tutorial to analyze multimodal using Seurat data step by step.

<!--more-->

### Metarial and Methods

- **Dataset**: 8617 cord blood mononuclear cells (CBMCs)


### Load in the data

- Valid link for the [RNA UMI matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz) download.

- Valid link for the [ADT UMI matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz) download.

    ```r
    library(Seurat)
    library(ggplot2)
    library(patchwork)

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

### Setup a Seurat object, add the RNA and protein data

- Create a Seurat object first, and then add the ADT data as a second assay

    ```r
    # creates a Seurat object based on the scRNA-seq data
    cbmc <- CreateSeuratObject(counts = cbmc.rna)

    # We can see that by default, the cbmc object contains an assay storing RNA measurement
    Assays(cbmc)

    # create a new assay to store ADT information
    adt_assay <- CreateAssayObject(counts = cbmc.adt)

    # add this assay to the previously created Seurat object
    cbmc[["ADT"]] <- adt_assay

    # Validate that the object now contains multiple assays
    Assays(cbmc)

    # Extract a list of features measured in the ADT assay
    rownames(cbmc[["ADT"]])

    # Note that we can easily switch back and forth between the two assays to specify the default
    # for visualization and analysis

    # List the current default assay
    DefaultAssay(cbmc)

    # Switch the default to ADT
    DefaultAssay(cbmc) <- "ADT"
    DefaultAssay(cbmc)
    ```

    ![fig2](2.png)

### Cluster cells on the basis of their scRNA-seq profiles

- The steps below represent a quick clustering of the PBMCs based on the scRNA-seq data. 

    ```r
    # Note that all operations below are performed on the RNA assay Set and verify that the default
    # assay is RNA
    DefaultAssay(cbmc) <- "RNA"
    DefaultAssay(cbmc)

    # perform visualization and clustering steps
    cbmc <- NormalizeData(cbmc)
    cbmc <- FindVariableFeatures(cbmc)
    cbmc <- ScaleData(cbmc)
    cbmc <- RunPCA(cbmc, verbose = FALSE)
    cbmc <- FindNeighbors(cbmc, dims = 1:30)
    cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
    cbmc <- RunUMAP(cbmc, dims = 1:30)
    DimPlot(cbmc, label = TRUE)
    ggsave("./dimplot1.png")
    ```

- Output:
  ![fig3](fig3.png)

- Dimplot result
  ![dimplot1](dimplot1.png)

### Visualize multiple modalities side-by-side

- Visualize the expression of either protein or RNA molecules in the dataset, after obtained clusters from scRNA-seq profiles.

- For example with the B cell marker CD19 (both protein and RNA levels).

    ```r
    # Normalize ADT data,
    DefaultAssay(cbmc) <- "ADT"
    cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
    DefaultAssay(cbmc) <- "RNA"

    # Note that the following command is an alternative but returns the same result
    cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

    # Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
    # visualize one or the other
    DefaultAssay(cbmc) <- "ADT"
    p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
    DefaultAssay(cbmc) <- "RNA"
    p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")

    # place plots side-by-side
    p1 | p2
    ggsave("./featureplot1.png")
    ```

    ![featureplot1](featureplot1.png)

- Herein, the merge and arrangement of generated plot using the `patchwork` package.

### Identify cell surface markers for scRNA-seq clusters


- Identify cluster 6 as expressing CD19 on the surface
  
    ```r
    VlnPlot(cbmc, "adt_CD19")
    ggsave("./vlnplot1.png")
    ```

    ![vlnplot1](vlnplot1.png)

- Also, identify alternative protein and RNA markers for this cluster through differential expression

    ```r
    adt_markers <- FindMarkers(cbmc, ident.1 = 5, assay = "ADT")
    rna_markers <- FindMarkers(cbmc, ident.1 = 5, assay = "RNA")

    head(adt_markers)
    head(rna_markers)
    ```

    ![fig3](fig3.png)

### Additional visualizations of multimodal data

- Draw `ADT` scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if desired by using HoverLocator and FeatureLocator.

    ```r  
    FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
    ggsave("./featurescatterplot.png")
    ```

    ![featurescatterplot](featurescatterplot.png)

- View relationship between **protein** and **RNA**.

    ```r
    FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
    ggsave("./featurescatterplot2.png")
    ```

    ![featurescatterplot2](featurescatterplot2.png)

- View relationship between **CD4 and CD8 proteins**.

    ```r
    FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")
    ggsave("./featurescatterplot3.png")
    ```

    ![featurescatterplot3](featurescatterplot3.png)

- View the **raw** (non-normalized) ADT `counts`.

    ```r
    FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")
    ggsave("./featurescatterplot4.png")
    ```

    ![featurescatterplot4](featurescatterplot4.png)

### Loading data from 10X multi-modal experiments

- **Download Dataset**: a dataset of 7,900 peripheral blood mononuclear cells, freely available from 10X Genomics [here](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.tar.gz).

- uncompress dataset.

    ```python
    tar xvzf pbmc_10k_protein_v3_filtered_feature_bc_matrix.tar.gz
    ```

- Analyze data.

    ```r
    pbmc10k.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
    rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", 
                                                            x = rownames(x = pbmc10k.data[["Antibody Capture"]]))

    pbmc10k <- CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 3, min.features = 200)
    pbmc10k <- NormalizeData(pbmc10k)
    pbmc10k[["ADT"]] <- CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
    pbmc10k <- NormalizeData(pbmc10k, assay = "ADT", normalization.method = "CLR")

    plot1 <- FeatureScatter(pbmc10k, feature1 = "adt_CD19", feature2 = "adt_CD3", pt.size = 1)
    plot2 <- FeatureScatter(pbmc10k, feature1 = "adt_CD4", feature2 = "adt_CD8a", pt.size = 1)
    plot3 <- FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "CD3E", pt.size = 1)
    (plot1 + plot2 + plot3) & NoLegend()
    ggsave("./featurescatter53.png")
    ```

    ![featurescatter53](featurescatter53.png)

### In summary

- There are additional functionality for multimodal data in Seurat. Keep exploring the resources.


