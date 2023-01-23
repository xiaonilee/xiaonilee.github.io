rm(list = ls())
options(stringsAsFactors = F)
# QC
## remove duplicated genes and do order

dat <- read.table("GSM3454528_naive_cells.txt.gz",
                   header = T,#row.names = 1,
                   stringsAsFactors = F)
head(dat)

dat <- dat[order(apply(dat[,-1], 1, sum),decreasing = T),]
dat <- dat[!duplicated(dat$genes),]
dat <- dat[order(dat$genes),]
rownames(dat) <- dat[,1]
dat <- dat[,-1]



## normalized

lib.size=colSums(dat)/median(colSums(dat))
dat.new=dat
for(i in 1:length(lib.size)){
  #print(i)
  dat.new[,i]=dat[,i]/lib.size[i]
}
dat=as.matrix(dat.new)


## log2()

dat=log2(dat+1)
dim(dat)


# choose top5000
## create Seurat

#install.packages(Seurat)
library("Seurat")
scRNA = CreateSeuratObject(counts=dat)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 5000) 
hvg.gene=VariableFeatures(scRNA)
str(hvg.gene)



dat.hvg=dat[rownames(dat) %in% hvg.gene,]
dim(dat.hvg)


## Do PCA with prcomp()

pca <-prcomp(t(dat.hvg))
dim(pca$x)
pca$x[1:4,1:4]



# k-means cluster
## k=10

clust.kmeans <- kmeans(pca$x[,1:20], centers=10)
table(clust.kmeans$cluster)


# KNN visualized

dist<-as.matrix(dist(pca$x[,1:20]))
dist[1:3,1:3]

edges <- mat.or.vec(0,2)


for (i in 1:nrow(dist)){
  # find closes neighbours
  matches <- setdiff(order(dist[i,],decreasing = F)[1:21],i) 
  
  # add edges
  edges <- rbind(edges,cbind(i,matches))  
}
head(edges, 50)


## Create igraph
 
library(igraph)
?igraph
graph <- graph_from_edgelist(edges,directed=F)
graph

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
 

## Define function(D,k) to enable to try different k
   
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
  #取消节点的边框颜色
  V(graph)$frame.color <- NA
  # make a layout for visualizing in 2D
  set.seed(1)
  #指定layout_with_fr类型布局风格
  g.layout<-layout_with_fr(graph)
  return(list(graph=graph,layout=g.layout))        
}
 


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
