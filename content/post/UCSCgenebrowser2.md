---
title: "UCSC gene browser with R package Gviz - part 2"
date: 2021-02-14
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

Continued, in this present, to draw the customized UCSC gene Track With Gviz package. 

<!--more-->

## Setup interested genomic location

```r
from <- 65921878
to <- 65980988
```

## Extract known genes in this region

```r
knownGenes <- UcscTrack(genome = "mm9", chromosome = "chrX",
                        track = "knownGene", from = from, to = to,
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds",
                        gene = "name", symbol = "name",
                        transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
```

## Extract reference genes from Ensemble database

```r
refGenes <- UcscTrack(genome = "mm9", chromosome = "chrX",
                      track = "xenoRefGene", from = from, to = to,
                      trackType = "GeneRegionTrack",
                      rstarts = "exonStarts", rends = "exonEnds",
                      gene = "name",  symbol = "name2",
                      transcript = "name", strand = "strand",
                      fill = "#8282d2", stacking = "dense",
                      name = "Other RefSeq")
ensGenes <- UcscTrack(genome = "mm9", chromosome = "chrX",
                      track = "ensGene", from = from, to = to,
                      trackType = "GeneRegionTrack",
                      rstarts = "exonStarts", rends = "exonEnds",
                      gene = "name", symbol = "name2",
                      transcript = "name", strand = "strand",
                      fill = "#960000", name = "Ensembl Genes")
```

## Extract location for CpGIsland and SNP

```r
cpgIslands <- UcscTrack(genome = "mm9", chromosome = "chrX",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack",
                        start = "chromStart", end = "chromEnd",
                        id = "name", shape = "box", fill = "#006400",
                        name = "CpG Islands")
snpLocations <-  UcscTrack(genome = "mm9", chromosome = "chrX",
                           track = "snp128", from = from, to = to,
                           trackType = "AnnotationTrack",
                           start = "chromStart", end = "chromEnd",
                           id = "name", feature = "func",
                           strand = "strand", shape = "box",
                           stacking = "dense", fill = "black",
                           name = "SNPs")
```

## conservation and GC content

```r
conservation <- UcscTrack(genome = "mm9", chromosome = "chrX",
                          track = "Conservation",
                          table = "phyloP30wayPlacental",
                          from = from, to = to, trackType = "DataTrack",
                          start = "start", end = "end", data = "score",
                          type = "hist", window = "auto",
                          col.histogram = "darkblue",
                          fill.histogram = "darkblue",
                          ylim = c(-3.7, 4), name = "Conservation")
gcContent <- UcscTrack(genome = "mm9", chromosome = "chrX",
                       track = "GC Percent", table = "gc5Base",
                       from = from, to = to, trackType = "DataTrack",
                       start = "start", end = "end", data = "score",
                       type = "hist", window = -1, windowSize = 1500,
                       fill.histogram = "black", col.histogram = "black",
                       ylim = c(30, 70), name = "GC Percent")
```

## lable information for aix and chromosome

```r
axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome="mm9", chromosome="chrX")
```

## Finally

```r
plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, ensGenes,
                cpgIslands, gcContent, conservation, snpLocations),
           from = from, to = to, showTitle = FALSE)
```

output
![fig1](fig1.png)
