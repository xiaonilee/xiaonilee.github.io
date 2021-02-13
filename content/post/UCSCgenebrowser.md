---
title: "UCSC gene browser with R package Gviz"
date: 2021-02-13
lastmod: 2021-02-13
draft: false
tags: ["R", "Bioinformatics", "UCSC"]
categories: ["R", "Bioinformatics", "UCSC"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Gviz, Visualize genomic data. Let's dive into it.

<!--more-->

## Introduction to Gviz package

- Provide a structured visualization framework to plot any type of data along genomic coordinates.
- Integrate publicly available genomic annotation data from sources like UCSC or ENSEMBL.

### Plot annotation track

```r
library(Gviz)
library(GenomicRanges)
#Load data : class = GRanges
data(cpgIslands)
cpgIslands
```

output: 
![p0](p0.png)

```r
#Annotation track, title ="CpG"
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)
```

output: 
![p1](p1.png)

### Add genome axis track

```r
## genomic coordinates
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))
```

output: 
![p2](p2.png)

### Add chromosome ideogram

```r
#genome : "hg19" 
gen<-genome(cpgIslands)
#Chromosme name : "chr7"
chr <- as.character(unique(seqnames(cpgIslands)))
#Ideogram track
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack))
```

output: 
![p3](p3.png)

### Add gene model

```r
#Load data
data(geneModels)
head(geneModels)
#Plot
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))
```

output: 
![p4](p4.png)

### Zoom the plot

```r
#Use from and to arguments to zoom
plotTracks(list(itrack, gtrack, atrack, grtrack),
           from = 26700000, to = 26750000)
```

output: 
![p5](p5.png)

```r
# Use extend.left and extend.right to zoom
#those arguments are relative to the currently displayed ranges, 
#and can be used to quickly extend the view on one or both ends of the plot.
plotTracks(list(itrack, gtrack, atrack, grtrack),
           extend.left = 0.5, extend.right = 1000000)
```

output: 
![p6](p6.png)

```r
# to drop the bounding borders of the exons and 
# to have a nice plot
plotTracks(list(itrack, gtrack, atrack, grtrack),
           extend.left = 0.5, extend.right = 1000000, col = NULL)
```

output:
![p7](p7.png)

### Add sequence track and zoom to view sequence

```r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack, grtrack,
                strack), from = 26591822, to = 26591852, cex = 0.8)
```

output:
![p8](p8.png)

### Add data track

```r
#For demonstration purposes we can create a simple DataTrack object from
#randomly sampled data.
set.seed(255)
lim <- c(26700000, 26750000)
coords <- sort(c(lim[1], sample(seq(from = lim[1], 
                                    to = lim[2]), 99), lim[2]))
dat <- runif(100, min = -10, max = 10)
head(dat)

##data track
dtrack <- DataTrack(data = dat, start = coords[-length(coords)],
                    end = coords[-1], chromosome = chr, genome = gen,
                    name = "Uniform")
##Plot data track
plotTracks(list(itrack, gtrack, atrack, grtrack,
                dtrack), from = lim[1], to = lim[2])
```

output:
![p9](p9.png)

```r
#Change plot type to histogram
plotTracks(list(itrack, gtrack, atrack, grtrack,dtrack),
           from = lim[1], to = lim[2], type = "histogram")
```

output:
![p10](p10.png)

### Plotting parameters

```r
#setting parameters
#Annotation of transcript
#Change panel and title background color
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")
plotTracks(list(itrack, gtrack, atrack, grtrack))
```

output:
![p11](p11.png)

### Plotting direction

```r
#By default all tracks will be plotted in a 5’ -> 3’ direction. 
plotTracks(list(itrack, gtrack, atrack, grtrack),
           reverseStrand = TRUE)
```

![p12](p12.png)

### Track classes

```r
#Set the position of labels to below, show IDs, change color


axisTrack <- GenomeAxisTrack(range = IRanges(start = c(2000000,4000000), 
                                             end = c(3000000, 7000000),
                                             names = rep("N-stretch", 2))
                             )
plotTracks(axisTrack, from = 1000000, to = 9000000, 
           labelPos = "below",showId=TRUE, col="red")
```

output:
![p13](p13.png)

### IdeogramTrack

```r
#Ideogram
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")
plotTracks(ideoTrack, from = 85000000, to = 129000000)
```

output:
![p14](p14.png)

```r
#Show chromosome band ID
plotTracks(ideoTrack, from = 85000000, to = 129000000,
           showId = FALSE, showBandId = TRUE, cex.bands = 0.4)
```

output:
![p15](p15.png)


### DataTrack

```r
#Load data
data(twoGroups)
head(twoGroups)
```

output
![p160](p160.png)

```r
#Plot data track
dTrack <- DataTrack(twoGroups, name = "uniform")
plotTracks(dTrack)
```

output:
![p16](p16.png)

```r
#The different plot types
#dotplot
plotTracks(DataTrack(twoGroups, name = "p"), type="p")
#lines plot
plotTracks(DataTrack(twoGroups, name = "l"), type="l")
#line and dot plot
plotTracks(DataTrack(twoGroups, name = "b"), type="b")
#lines plot of average
plotTracks(DataTrack(twoGroups, name = "a"), type="a")
#histogram lines
plotTracks(DataTrack(twoGroups, name = "h"), type="h")
#histogram histogram (bar width equal to range with)
plotTracks(DataTrack(twoGroups, name = "histogram"), type="histogram")
#'polygon-type' plot relative to a baseline
plotTracks(DataTrack(twoGroups, name = "polygon"), type="polygon")
#box and whisker plot
plotTracks(DataTrack(twoGroups, name = "boxplot"), type="boxplot")
#false color image of the individual values
plotTracks(DataTrack(twoGroups, name = "heatmap"), type="heatmap")
```

output
![p17](p17.png)

### Example of DataTrack plots

```r
#Combine a boxplot with an average line and a data grid (g):
plotTracks(dTrack, type = c("boxplot", "a", "g"))
#Heatmap and show sample names
plotTracks(dTrack, type = c("heatmap"), showSampleNames = TRUE,
           cex.sampleNames = 0.6)
```

output
![p18](p18.png)

### Data grouping

```r
plotTracks(dTrack, groups = rep(c("control", "treated"), each = 3),
           type = c("a", "p"), legend=TRUE)
#Boxplot
plotTracks(dTrack, groups = rep(c("control", "treated"), each = 3),
           type = "boxplot")
#Aggregate group. aggregation can be  mean, median, extreme,
#sum, min and max
plotTracks(dTrack, groups = rep(c("control", "treated"),each = 3),
           type = c("b"), aggregateGroups = TRUE,
           aggregation = "max")
```
