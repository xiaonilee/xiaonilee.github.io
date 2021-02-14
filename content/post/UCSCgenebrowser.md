---
title: "UCSC gene browser with R package Gviz - part 1"
date: 2021-02-13
lastmod: 2021-02-14
draft: false
tags: ["R", "Bioinformatics", "UCSC", "R Package"]
categories: ["R", "Bioinformatics", "UCSC", "R Package"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Gviz, Visualize genomic data. Let's dive into it to explore useful functions.

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

output 
![p0](p0.png)

```r
#Annotation track, title ="CpG"
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)
```

output 
![p1](p1.png)

### Add genome axis track

```r
## genomic coordinates
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))
```

output 
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

output 
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

output 
![p4](p4.png)

### Zoom the plot

```r
#Use from and to arguments to zoom
plotTracks(list(itrack, gtrack, atrack, grtrack),
           from = 26700000, to = 26750000)
```

output 
![p5](p5.png)

```r
# Use extend.left and extend.right to zoom
#those arguments are relative to the currently displayed ranges, 
#and can be used to quickly extend the view on one or both ends of the plot.
plotTracks(list(itrack, gtrack, atrack, grtrack),
           extend.left = 0.5, extend.right = 1000000)
```

output 
![p6](p6.png)

```r
# to drop the bounding borders of the exons and 
# to have a nice plot
plotTracks(list(itrack, gtrack, atrack, grtrack),
           extend.left = 0.5, extend.right = 1000000, col = NULL)
```

output
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

output
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

output
![p9](p9.png)

```r
#Change plot type to histogram
plotTracks(list(itrack, gtrack, atrack, grtrack,dtrack),
           from = lim[1], to = lim[2], type = "histogram")
```

output
![p10](p10.png)

## Plotting parameters

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

output
![p11](p11.png)

### Plotting direction

```r
#By default all tracks will be plotted in a 5’ -> 3’ direction. 
plotTracks(list(itrack, gtrack, atrack, grtrack),
           reverseStrand = TRUE)
```

![p12](p12.png)

## Track classes

```r
#Set the position of labels to below, show IDs, change color


axisTrack <- GenomeAxisTrack(range = IRanges(start = c(2000000,4000000), 
                                             end = c(3000000, 7000000),
                                             names = rep("N-stretch", 2))
                             )
plotTracks(axisTrack, from = 1000000, to = 9000000, 
           labelPos = "below",showId=TRUE, col="red")
```

output
![p13](p13.png)

### IdeogramTrack

```r
#Ideogram
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")
plotTracks(ideoTrack, from = 85000000, to = 129000000)
```

output
![p14](p14.png)

```r
#Show chromosome band ID
plotTracks(ideoTrack, from = 85000000, to = 129000000,
           showId = FALSE, showBandId = TRUE, cex.bands = 0.4)
```

output
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

output
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
![p17](p17new.png)

### Example of DataTrack plots

```r
#Combine a boxplot with an average line and a data grid (g):
plotTracks(dTrack, type = c("boxplot", "a", "g"))
#Heatmap and show sample names
plotTracks(dTrack, type = c("heatmap"), showSampleNames = TRUE,
           cex.sampleNames = 0.6)
```

output
![p18](p18new.png)

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

output
![p19](p19.png)

### Building DataTrack objects from files

```r
bgFile <- system.file("extdata/test.bedGraph", package = "Gviz")
dTrack2 <- DataTrack(range = bgFile, genome = "hg19",
                     type = "l", chromosome = "chr19", name = "bedGraph")
plotTracks(dTrack2)

bamFile <- system.file("extdata/test.bam", package = "Gviz")
dTrack4 <- DataTrack(range = bamFile, genome = "hg19",
                     type = "l", name = "Coverage", window = -1, chromosome = "chr1")
plotTracks(dTrack4, from = 189990000, to = 190000000)
```

output
![p20](p20new.png)

### AnnotationTrack

```r
aTrack <- AnnotationTrack(start = c(10, 40, 120),
                          width = 15, chromosome = "chrX", strand = c("+","*", "-"),
                          id = c("Huey", "Dewey", "Louie"),
                          genome = "hg19", name = "foo")
plotTracks(aTrack)
```

output
![p21](p21.png)

### Add id information

```r
plotTracks(aTrack, shape = "box", featureAnnotation = "id")
```

output
![p210](p210.png)

### change fontcolor

```r
plotTracks(aTrack, shape = "ellipse", featureAnnotation = "id", 
          fontcolor.feature = "darkblue")
### Building AnnotationTrack objects from files
```

output
![p211](p211.png)



```r
#Annotation track
aTrack2 <- AnnotationTrack(range = bamFile, genome = "hg19",
                           name = "Reads", chromosome = "chr1")
plotTracks(aTrack2, from = 189995000, to = 190000000)
```

output
![p22](p22.png)


#### plot both the DataTrack representation as well as the AnnotationTrack representation of the bam file together

```r

plotTracks(list(dTrack4, aTrack2), from = 189990000,
           to = 190000000)
```

output
![p23](p23.png)

### GeneRegionTrack

```r
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "foo", 
                           transcriptAnnotation = "symbol")

```

### Building GeneRegionTrack objects from TranscriptDbs

```r
library(GenomicFeatures)
samplefile <- system.file("extdata/", "UCSC_knownGene_sample.sqlite",
                          package = "GenomicFeatures")
txdb <- loadDb(samplefile,packageName = "GenomicFeatures")
txTr <- GeneRegionTrack(txdb, chromosome = "chr6", start = 300000, end = 350000)
#feature(txTr)
plotTracks(txTr)
```

### BiomartGeneRegionTrack

```r
biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                    chromosome = chr, start = 20000000, end =
                                      21000000,name = "ENSEMBL")
plotTracks(biomTrack, col.line = NULL, col = NULL)
```

output
![p24](p24.png)


### Sequence Track

```r
library(BSgenome.Hsapiens.UCSC.hg19)
sTrack <- SequenceTrack(Hsapiens)
#sequence track : add 5'->3'
plotTracks(sTrack, chromosome = 1, from = 20000,to = 20050,
           add53=TRUE)
#The complement
plotTracks(sTrack, chromosome = 1, from = 20000,to = 20050,
           add53=TRUE, complement = TRUE)
```

output
![p25](p25.png)

### AlignmentsTrack
#### RNAseq experiment

```r
afrom=2960000
ato=3160000
#bam file
alTrack <- AlignmentsTrack(system.file(package = "Gviz", "extdata", "gapped.bam"), isPaired = TRUE)
bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr12",
                              start = afrom, end = ato, filter = list(with_ox_refseq_mrna = TRUE),
                              stacking = "dense")
plotTracks(c(bmt, alTrack), from = afrom, to = ato, chromosome = "chr12")
```

output
![p26](p26.png)


- To reduce the size of the coverage section by setting the coverageHeight or the minCoverageHeight parameters.

```r
plotTracks(c(bmt, alTrack), from = afrom, to = ato,
           chromosome = "chr12", min.height = 0, coverageHeight = 0.08,
           minCoverageHeight = 0)
```

output
![p27](p27.png)


- Turn off the pile-up view of the individual reads by setting the type display parameter.

```r
plotTracks(c(alTrack, bmt), from = afrom, to = ato, chromosome = "chr12", type = "coverage")
```

output
![p28](p28.png)

### Zoom in a bit further to check out the details of the pile-ups section

```r
plotTracks(c(bmt, alTrack), from = afrom + 12700,
           to = afrom + 15200, chromosome = "chr12")
```

output
![p29](p29.png)


- in single end mode by setting the isPaired argument in the constructor.

```r
alTrack <- AlignmentsTrack(system.file(package = "Gviz",
                                       "extdata", "gapped.bam"), isPaired = FALSE)
plotTracks(c(bmt, alTrack), from = afrom + 12700,
           to = afrom + 15200, chromosome = "chr12")
```

output
![p30](p30.png)

### DNAseq experiment

```r
afrom <- 44945200
ato <- 44947200
alTrack <- AlignmentsTrack(system.file(package = "Gviz","extdata", "snps.bam"), isPaired = TRUE)
plotTracks(c(alTrack, sTrack), chromosome = "chr21", from = afrom,to = ato)
```

output
![p31](p31.png)

- Zoom in to one of the obvious heterozygous SNP positions 

```r
plotTracks(c(alTrack, sTrack), chromosome = "chr21", from = 44946590, to = 44946660)
#show individual letters
plotTracks(c(alTrack, sTrack), chromosome = "chr21",
           from = 44946590, to = 44946660, cex = 0.5, min.height = 8)
```

output
![p32](p32.png)

## Track highlighting and overlays

```r
#highlight
ht <- HighlightTrack(trackList = list(atrack, grtrack,dtrack),
                     start = c(26705000, 26720000), width = 7000, chromosome = 7)
plotTracks(list(itrack, gtrack, ht), from = lim[1], to = lim[2])
```

output
![p33](p33.png)

### overlays multiple tracks on the same area of the plot

```r
#create data
dat <- runif(100, min = -2, max = 22)
dtrack2 <- DataTrack(data = dat, start = coords[-length(coords)],
                     end = coords[-1], chromosome = chr, genome = gen,
                     name = "Uniform2", groups = factor("sample 2",levels = c("sample 1", "sample 2")),
                     legend = TRUE)
displayPars(dtrack) <- list(groups = factor("sample 1",levels = c("sample 1", "sample 2")), legend = TRUE)
ot <- OverlayTrack(trackList = list(dtrack2, dtrack))
ylims <- extendrange(range(c(values(dtrack), values(dtrack2))))
plotTracks(list(itrack, gtrack, ot), from = lim[1], to = lim[2], ylim = ylims, type = c("smooth", "p"))
```

output
![p34](p34.png)


- alpha blending can be a useful tool to tease out even more information out of track overlays

```r
displayPars(dtrack) <- list(alpha.title = 1, alpha = 0.5)
displayPars(dtrack2) <- list(alpha.title = 1, alpha = 0.5)
ot <- OverlayTrack(trackList = list(dtrack, dtrack2))
plotTracks(list(itrack, gtrack, ot), from = lim[1],
           to = lim[2], ylim = ylims, type = c("hist"), window = 30)
```

output
![p35](p35.png)
