
# Chapter 16. Cluster analysis


# This chapter covers
#   Identifying cohesive subgroups (clusters) of observations
#   Determining the number of clusters present
#   Obtaining a nested hierarchy of clusters
#   Obtaining discrete clusters


# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 16.2. Calculating distances
# install.packages("flexclust")
data(nutrient, package = "flexclust")
head(nutrient, 4)
#              energy protein fat calcium iron
# BEEF BRAISED    340      20  28       9  2.6
# HAMBURGER       245      21  17       9  2.7
# BEEF ROAST      420      15  39       7  2.0
# BEEF STEAK      375      19  32       9  2.6

# options(digits = 4)
d <- dist(nutrient)
as.matrix(d)[1:4, 1:4]
#              BEEF BRAISED HAMBURGER BEEF ROAST BEEF STEAK
# BEEF BRAISED         0.00     95.64      80.93      35.24
# HAMBURGER           95.64      0.00     176.49     130.88
# BEEF ROAST          80.93    176.49       0.00      45.76
# BEEF STEAK          35.24    130.88      45.76       0.00


# 16.3 Hierarchical cluster analysis
# code listing 16.1. Average-linkage clustering of the nutrient data
data(nutrient, package = "flexclust")
row.names(nutrient) <- tolower(row.names(nutrient))
nutrient.scaled <- scale(nutrient)

d <- dist(nutrient.scaled)
fit.average <- hclust(d, method = "average")
# figure 16.1
plot(fit.average, hang = -1, cex=0.8, main = "Average Linkage Clustering")  


#=================================================================================
# code listing 16.2. Selecting the number of clusters
# install.packages("NbClust")
library(NbClust)
devAskNewPage(ask = F)
nc <- NbClust(nutrient.scaled, distance = "euclidean",
              min.nc = 2, max.nc = 15, method = "average")

table(nc$Best.nc[1,])
# 0  1  2  3  4  5  9 10 13 14 15 
# 2  1  4  4  2  4  1  1  2  1  4

barplot(table(nc$Best.nc[1,]),
        xlab = "Number of Cluster", ylab = "Number of Criteria",
        main = "Number of Cluster Chosen by 26 Criteria") # figure 16.2

#=================================================================================
# code listing 16.3. Obtaining the final cluster solution
clusters <- cutree(fit.average, k=5)
table(clusters)
# clusters
# 1  2  3  4  5 
# 7 16  1  2  1

aggregate(nutrient, by=list(cluster=clusters), median)
#   cluster energy protein fat calcium iron
# 1       1  340.0      19  29       9 2.50
# 2       2  170.0      20   8      13 1.45
# 3       3  160.0      26   5      14 5.90
# 4       4   57.5       9   1      78 5.70
# 5       5  180.0      22   9     367 2.50

aggregate(as.data.frame(nutrient.scaled), by=list(cluster=clusters), median)
#   cluster  energy protein     fat calcium    iron
# 1       1  1.3101  0.0000  1.3786 -0.4480  0.0811
# 2       2 -0.3696  0.2352 -0.4869 -0.3968 -0.6374
# 3       3 -0.4684  1.6464 -0.7534 -0.3840  2.4078
# 4       4 -1.4812 -2.3520 -1.1088  0.4362  2.2709
# 5       5 -0.2708  0.7056 -0.3981  4.1397  0.0811

plot(fit.average, hang = -1, cex=0.8,
     main = "Average Linkage Clustering\n5 Cluster Solution") # figure 16.3-1

rect.hclust(fit.average, k=5) # figure 16.3-2

# 16.4.1 K-means clustering
wssplot <- function(data, nc=15, seed=1234) {
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc) {
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss) }
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares") 
  }

# code listing 16.4. K-means clustering of wine data
# install.packages("rattle")
data(wine, package = "rattle")
head(wine)

#   Type Alcohol Malic  Ash Alcalinity Magnesium Phenols Flavanoids Nonflavanoids
# 1    1   14.23  1.71 2.43       15.6       127    2.80       3.06          0.28
# 2    1   13.20  1.78 2.14       11.2       100    2.65       2.76          0.26
# 3    1   13.16  2.36 2.67       18.6       101    2.80       3.24          0.30
# 4    1   14.37  1.95 2.50       16.8       113    3.85       3.49          0.24
# 5    1   13.24  2.59 2.87       21.0       118    2.80       2.69          0.39
# 6    1   14.20  1.76 2.45       15.2       112    3.27       3.39          0.34
#   Proanthocyanins Color  Hue Dilution Proline
# 1            2.29  5.64 1.04     3.92    1065
# 2            1.28  4.38 1.05     3.40    1050
# 3            2.81  5.68 1.03     3.17    1185
# 4            2.18  7.80 0.86     3.45    1480
# 5            1.82  4.32 1.04     2.93     735
# 6            1.97  6.75 1.05     2.85    1450

df <- scale(wine[-1])
wssplot(df)
library(NbClust)
set.seed(1234)
devAskNewPage(ask = F)
nc <- NbClust(df, min.nc = 2, max.nc = 15, method = "kmeans")
table(nc$Best.nc[1,])

barplot(table(nc$Best.nc[1,]),
        xlab="Number of Clusters", ylab="Number of Criteria",
        main="Number of Clusters Chosen by 26 Criteria")
set.seed(1234)
fit.km <- kmeans(df, 3, nstart = 25)
fit.km$size
# [1] 62 65 51

fit.km$centers
#      Alcohol      Malic        Ash Alcalinity   Magnesium     Phenols  Flavanoids
# 1  0.8328826 -0.3029551  0.3636801 -0.6084749  0.57596208  0.88274724  0.97506900
# 2 -0.9234669 -0.3929331 -0.4931257  0.1701220 -0.49032869 -0.07576891  0.02075402
# 3  0.1644436  0.8690954  0.1863726  0.5228924 -0.07526047 -0.97657548 -1.21182921
#   Nonflavanoids Proanthocyanins      Color        Hue   Dilution    Proline
# 1   -0.56050853      0.57865427  0.1705823  0.4726504  0.7770551  1.1220202
# 2   -0.03343924      0.05810161 -0.8993770  0.4605046  0.2700025 -0.7517257
# 3    0.72402116     -0.77751312  0.9388902 -1.1615122 -1.2887761 -0.4059428

aggregate(wine[-1], by=list(cluster=fit.km$cluster), mean)
#   cluster  Alcohol    Malic      Ash Alcalinity Magnesium  Phenols Flavanoids Nonflavanoids
# 1       1 13.67677 1.997903 2.466290   17.46290 107.96774 2.847581  3.0032258     0.2920968
# 2       2 12.25092 1.897385 2.231231   20.06308  92.73846 2.247692  2.0500000     0.3576923
# 3       3 13.13412 3.307255 2.417647   21.24118  98.66667 1.683922  0.8188235     0.4519608
#   Proanthocyanins    Color       Hue Dilution   Proline
# 1        1.922097 5.453548 1.0654839 3.163387 1100.2258
# 2        1.624154 2.973077 1.0627077 2.803385  510.1692
# 3        1.145882 7.234706 0.6919608 1.696667  619.0588

ct.km <- table(wine$Type, fit.km$cluster)
ct.km
#    1  2  3
# 1 59  0  0
# 2  3 65  3
# 3  0  0 48

library(flexclust)
randIndex(ct.km)
#      ARI 
# 0.897495 

# 16.4.2. Partitioning around medoids
# code listing 16.5. Partitioning around medoids for the wine data
library(cluster)
set.seed(1234)
fit.pam <- pam(wine[-1], k=3, stand = TRUE)
fit.pam$medoids
#      Alcohol Malic  Ash Alcalinity Magnesium Phenols Flavanoids Nonflavanoids
# [1,]   13.48  1.81 2.41       20.5       100    2.70       2.98          0.26
# [2,]   12.25  1.73 2.12       19.0        80    1.65       2.03          0.37
# [3,]   13.40  3.91 2.48       23.0       102    1.80       0.75          0.43
#      Proanthocyanins Color  Hue Dilution Proline
# [1,]            1.86   5.1 1.04     3.47     920
# [2,]            1.63   3.4 1.00     3.17     510
# [3,]            1.41   7.3 0.70     1.56     750

clusplot(fit.pam, main = "Bivariate Cluster Plot")  # figure 16.6

ct.pam <- table(wine$Type, fit.pam$clustering)
ct.pam
#    1  2  3
# 1 59  0  0
# 2 16 53  2
# 3  0  1 47
randIndex(ct.pam)
#       ARI 
# 0.6994957

# 16.5. Avoiding nonexistent clusters
# figure 16.7
# install.packages("fMultivar")
library(fMultivar)
set.seed(1234)
df <- rnorm2d(1000, rho = .5)
df <- as.data.frame(df)
plot(df, main = "Bivariate Normal Distribution with rho=0.5")

# figure 16.8
wssplot(df)
library(NbClust)
nc <- NbClust(df, min.nc = 2, max.nc = 15, method = "kmeans")
dev.new()

# figure 16.9
barplot(table(nc$Best.nc[1,]),
        xlab="Number of Clusters", ylab="Number of Criteria",
        main="Number of Clusters Chosen by 26 Criteria")


library(ggplot2)
library(cluster)
fit <- pam(df, k=2)
df$clustering <- factor(fit$clustering)
# figure 16.10
ggplot(data=df, aes(x=V1, y=V2, color=clustering, shape=clustering)) +
  geom_point() + ggtitle("Clustering of Bivariate Normal Data")

# figure 16.11
plot(nc$All.index[,4], type="o", ylab="CCC",
     xlab="Number of clusters", col="blue")