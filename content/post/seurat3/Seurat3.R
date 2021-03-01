rm(list = ls())
options(stringsAsFactors = F)

#================10x Visium=================

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


#========Load dataset==========
#BiocManager::install("rhdf5")
#BiocManager::install("rhdf5r")
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

#==========Data preprocessing=================
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
ggsave("./fig1.png")

#sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

#==================Gene expression visualization=======================

#SpatialFeaturePlot() function overlay molecular data on top of tissue histology. 
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
ggsave("./fig2.png")

p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2
ggsave("./fig3.png")


#============Dimensionality reduction, clustering, and visualization========

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
ggsave("./fig4.png")

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)

ggsave("./fig5.png")

#==============Interactive plotting==============
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)
ggsave("./fig6.png")

#=========Identification of Spatially Variable Features=======
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, 
                   features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
ggsave("./fig7.png")


# brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", 
#                                        features = VariableFeatures(brain)[1:1000], 
#                                        selection.method = "markvariogram")
# ggsave("./fig8.png")


# top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
# SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
# ggsave("./fig9.png")


#=============Subset out anatomical regions====================

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


#==========Integration with single-cell data==============
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

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE, 
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay


DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE)
ggsave("./fig12.png")


cortex <- FindSpatiallyVariableFeatures(cortex, 
                                        assay = "predictions", 
                                        selection.method = "markvariogram", 
                                        features = rownames(cortex), 
                                        r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
ggsave("./fig13.png")


SpatialFeaturePlot(cortex, 
                   features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", 
                                        "L6b", "Oligo"), 
                   pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
ggsave("./fig14.png")


#============Working with multiple slices in Seurat==================

brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

brain.merge <- merge(brain, brain2)

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)


DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
ggsave("./fig15.png")

SpatialDimPlot(brain.merge)
ggsave("./fig16.png")

SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))
ggsave("./fig17.png")

#===================Slide-seq=================

InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")

plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
ggsave("./fig18.png")

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)


plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2
ggsave("./fig19.png")

SpatialDimPlot(slide.seq, 
               cells.highlight = CellsByIdentities(object = slide.seq, 
                                                   idents = c(1, 6, 13)), 
               facet.highlight = TRUE)
ggsave("./fig20.png")                                            
                                                                                             