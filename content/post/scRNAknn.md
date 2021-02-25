---
title: "scRNA-sequencing analysis - PCA and KNN"
date: 2021-02-23
lastmod: 2021-02-24
draft: false
tags: ["Bioinformatics", "R", "RNA Sequence", "PCA", "KNN"]
categories: ["Bioinformatics", "R", "RNA Sequence", "PCA", "KNN"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---



<!--more-->

### Material and methods

- **Data set**: `GSE122083` / `GSM3454528`

- **R package**: `Seurat` and `igraph`

### Step 1. QC

- Download [data](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3454528&format=file&file=GSM3454528%5Fnaive%5Fcells%2Etxt%2Egz)

- Remove duplicated genes and Select max expression

```R
dat <- read.table("GSM3454528_naive_cells.txt.gz",
                   header = T,#row.names = 1,
                   stringsAsFactors = F)
head(dat)

dat <- dat[order(apply(dat[,-1], 1, sum),decreasing = T),]
dat <- dat[!duplicated(dat$genes),]
dat <- dat[order(dat$genes),]
rownames(dat) <- dat[,1]
dat <- dat[,-1]
```

- Normalized

```R
lib.size=colSums(dat)/median(colSums(dat))
dat.new=dat
for(i in 1:length(lib.size)){
  #print(i)
  dat.new[,i]=dat[,i]/lib.size[i]
}
dat=as.matrix(dat.new)
```

- Transform with log2()

```R
dat=log2(dat+1)
dim(dat)
```

  ![fig1](fig1.png)

### Step 2. Choose top5000 genes for PCA with `Seurat` package

- [Seurat](https://satijalab.org/seurat/index.html)
  ![fig8](fig8.png)

- Create Seurat

```R
#install.packages(Seurat)
library("Seurat")
scRNA = CreateSeuratObject(counts=dat)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 5000) 
hvg.gene=VariableFeatures(scRNA)
str(hvg.gene)
```

  ![fig2](fig2.png)

- PCA with `prcomp()`

```R
pca <-prcomp(t(dat.hvg))
dim(pca$x)
pca$x[1:4,1:4]
```

  ![fig3](fig3.png)

### Step 3. k-means cluster

```R
# k=10
clust.kmeans <- kmeans(pca$x[,1:20], centers=10)
table(clust.kmeans$cluster)
```

  ![fig4](fig4.png)

### Step 4. KNN visualized

```R
dist<-as.matrix(dist(pca$x[,1:20]))
dist[1:3,1:3]

edges <- mat.or.vec(0,2)
```

  ![fig5](fig5.png)

```R
# iter. for every cell
for (i in 1:nrow(dist)){
  # find closes neighbours
  matches <- setdiff(order(dist[i,],decreasing = F)[1:21],i) 
  
  # add edges
  edges <- rbind(edges,cbind(i,matches))  
}
head(edges, 50)
```

- Create igraph

```R
library(igraph)
?igraph
graph <- graph_from_edgelist(edges,directed=F)
graph
```

  ![fig6](fig6.png)

```R
# lable colors for cells.
cols<-rainbow(10)
names(cols) <- unique(clust.kmeans$cluster)
col.clust <- cols[clust.kmeans$cluster]

png("test1.png")
set.seed(1)
plot(graph,vertex.size=1,vertex.label=NA,vertex.frame.color=NA,vertex.color=col.clust,
            edge.width=0.5,main="50PCs; k=20")
legend("topright",names(cols),col=cols,
       pch=16,cex=0.5,bty='n')
dev.off()
```

  ![test1](test1.png)


- Define function(D,k) to enable to try different k

```R
make.knn.graph<-function(D,k){
  # calculate euclidean distances between cells
  dist<-as.matrix(dist(D))
  # make a list of edges to k nearest neighbors for each cell
  edges <- mat.or.vec(0,2)
  for (i in 1:nrow(dist)){
    # find closes neighbours
    matches <- setdiff(order(dist[i,],decreasing = F)[1:(k+1)],i)
    # add edges in both directions
    edges <- rbind(edges,cbind(rep(i,k),matches))  
  }
  # create a graph from the edgelist
  graph <- graph_from_edgelist(edges,directed=F)
  # set frame color to null
  V(graph)$frame.color <- NA
  # make a layout for visualizing in 2D
  set.seed(1)
  # set layout to layout_with_fr
  g.layout<-layout_with_fr(graph)
  return(list(graph=graph,layout=g.layout))        
}
```

- Plot with different k values

```R
ks <- c(5,10,30,50)
for (k in ks){
  g.pca20 <- make.knn.graph(pca$x[,1:20],k)
  # plot all 4:
  print(k)
  png(paste0("k",k,".png"))
  plot.igraph(g.pca20$graph,layout=g.pca20$layout,vertex.color=col.clust,
            vertex.size=1,vertex.label=NA,main=paste0("K",k,"--","50 PCs"))
  legend("topright",names(cols),col=cols,
       pch=16,cex=0.5,bty='n')
  dev.off()
}
```

  ![fig7](fig7.png)


### In summary

- PCA and KNN are common methods, herein, they are standard ways to implement.

- Keep different parameters could get different results.

