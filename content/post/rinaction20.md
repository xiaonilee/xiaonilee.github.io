---
title: "Chapter 20. Advanced programming"
date: 2020-11-19
lastmod: 2020-11-19
draft: false
tags: ["R", "R in Action", "Bioinformatics", "Book"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

Notebook of Reading Books: R in Action_Chapter 20.

<!--more-->

## This chapter covers

- A deeper dive into the R language

- Using R’s OOP features to create generic functions

- Tweaking code to run more efficiently

- Finding and correcting programming errors

### 20.1. A review of the language

#### 20.1.1. Data types

- Atomic Vectors
- Generic Vectors or Lists
- Indexing

- Figure 20.1. A plot of the centroids (means) for three clusters extracted from the Iris dataset using ***k-means clustering***.
  
  ![fig201](fig201.png)

#### 20.1.2. Control structures

- `for()`
- `if()`
- `ifelse()`

#### 20.1.3. Creating functions

- Function syntax
- Object scope

### 20.2. Working with environments

### 20.3. Object-oriented programming

### 20.4. Writing efficient code

- Efficient data input
- Vectorization
- Correctly sizing objects
- Parallelization

### 20.5. Debugging

#### 20.5.1. Common sources of errors

#### 20.5.2. Debugging tools

- Built-in debugging functions

  ![tab1](tab1.png)

#### 20.5.3. Session options that support debugging

- options(error=traceback)
- options(error=recover)

### 20.6. Going further

There are a number of excellent sources of information on advanced programming in R.

- "**The R Language Definition**" (<http://mng.bz/U4Cm>) is a good place to ***start***. 

- “**Frames, Environments, and Scope in R and S-PLUS**” by John Fox (<http://mng.bz/Kkbi>) is a great article for gaining a better understanding of ***scope***. 

- “**How R Searches and Finds Stuff**” by Suraj Gupta (<http://mng.bz/2o5B>), is a blog article that can help you understand just what the title implies(`Going further`). 

- To learn more about ***efficient coding***, see “**FasteR! HigheR! StrongeR!—A Guide to Speeding Up R Code for Busy People**” by Noam Ross (<http://mng.bz/Iq3i>). 

- Finally, **R Programming for Bioinformatics** (<http://cran.r-project.org/web/views/ExperimentalDesign.html.>) by Robert Gentleman is a comprehensive yet highly accessible text for programmers that want to look under the hood. 

Attach is the [Script](chapter20.R) of chapter20.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 20.1. A review of the language

# 20.1.1. Data types

## Atomic Vectors
passed <- c(TRUE, TRUE, FALSE, TRUE)
passed
# [1]  TRUE  TRUE FALSE  TRUE

ages <- c(15, 18, 25, 14, 19)
ages
# [1] 15 18 25 14 19

cmplxNums <- c(1+2i, 0+1i, 39+3i, 12+2i)
cmplxNums
# [1]  1+2i  0+1i 39+3i 12+2i

names <- c("Bob", "Ted", "Carol", "Alice")
names
# [1] "Bob"   "Ted"   "Carol" "Alice"

x <- c(1,2,3,4,5,6,7,8)
class(x)
# [1] "numeric"

x
# [1] 1 2 3 4 5 6 7 8

attr(x, "dim") <- c(2, 4)
x
#      [,1] [,2] [,3] [,4]
# [1,]    1    3    5    7
# [2,]    2    4    6    8

class(x)
# [1] "matrix" "array" 

attributes(x)
# $dim
# [1] 2 4

attr(x, "dimnames") <- list(c("A1", "A2"),
                            c("B1", "B2", "B3", "B4"))
x
#    B1 B2 B3 B4
# A1  1  3  5  7
# A2  2  4  6  8

attr(x, "dim") <- NULL
class(x)
# [1] "numeric"

x
# [1] 1 2 3 4 5 6 7 8

## Generic Vectors or Lists
head(iris)
#   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
# 1          5.1         3.5          1.4         0.2  setosa
# 2          4.9         3.0          1.4         0.2  setosa
# 3          4.7         3.2          1.3         0.2  setosa
# 4          4.6         3.1          1.5         0.2  setosa
# 5          5.0         3.6          1.4         0.2  setosa
# 6          5.4         3.9          1.7         0.4  setosa

unclass(iris)
attributes(iris)
str(iris)

set.seed(1234)
fit <- kmeans(iris[1:4], 3)
names(fit)
# [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
# [6] "betweenss"    "size"         "iter"         "ifault"


unclass(fit)
# $cluster
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# [44] 1 1 1 1 1 1 1 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2
# [87] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 3 3 3 3 2 3 3 3 3 3 3 2 2 3 3 3 3 2 3 2 3 2 3 3 2 2 3
# [130] 3 3 3 3 2 3 3 3 3 2 3 3 3 2 3 3 3 2 3 3 2
# 
# $centers
# Sepal.Length Sepal.Width Petal.Length Petal.Width
# 1     5.006000    3.428000     1.462000    0.246000
# 2     5.901613    2.748387     4.393548    1.433871
# 3     6.850000    3.073684     5.742105    2.071053
# 
# $totss
# [1] 681.3706
# 
# $withinss
# [1] 15.15100 39.82097 23.87947
# 
# $tot.withinss
# [1] 78.85144
# 
# $betweenss
# [1] 602.5192
# 
# $size
# [1] 50 62 38
# 
# $iter
# [1] 2
# 
# $ifault
# [1] 0

sapply(fit, class)
# $cluster
# [1] "integer"
# 
# $centers
# [1] "matrix" "array" 
# 
# $totss
# [1] "numeric"
# 
# $withinss
# [1] "numeric"
# 
# $tot.withinss
# [1] "numeric"
# 
# $betweenss
# [1] "numeric"
# 
# $size
# [1] "integer"
# 
# $iter
# [1] "integer"
# 
# $ifault
# [1] "integer"


#===========================================
## indexing
x <- c(20, 30, 40)
x[3]
# [1] 40

x[c(2, 3)]
# [1] 30 40

x <- c(A=20, B=30, C=40)
x
#  A  B  C 
# 20 30 40

x[c(2,3)]
#  B  C 
# 30 40 

x[c("B", "C")]
#  B  C 
# 30 40 

str(fit)
# List of 9
# $ cluster     : int [1:150] 1 1 1 1 1 1 1 1 1 1 ...
# $ centers     : num [1:3, 1:4] 5.01 5.9 6.85 3.43 2.75 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "1" "2" "3"
# .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
# $ totss       : num 681
# $ withinss    : num [1:3] 15.2 39.8 23.9
# $ tot.withinss: num 78.9
# $ betweenss   : num 603
# $ size        : int [1:3] 50 62 38
# $ iter        : int 2
# $ ifault      : int 0
# - attr(*, "class")= chr "kmeans"

fit[7]
# $size
# [1] 50 62 38


fit[c(2,7)]
# $centers
#   Sepal.Length Sepal.Width Petal.Length Petal.Width
# 1     5.006000    3.428000     1.462000    0.246000
# 2     5.901613    2.748387     4.393548    1.433871
# 3     6.850000    3.073684     5.742105    2.071053
# 
# $size
# [1] 50 62 38


fit[2]
# $centers
#   Sepal.Length Sepal.Width Petal.Length Petal.Width
# 1     5.006000    3.428000     1.462000    0.246000
# 2     5.901613    2.748387     4.393548    1.433871
# 3     6.850000    3.073684     5.742105    2.071053

fit[[2]]
#   Sepal.Length Sepal.Width Petal.Length Petal.Width
# 1     5.006000    3.428000     1.462000    0.246000
# 2     5.901613    2.748387     4.393548    1.433871
# 3     6.850000    3.073684     5.742105    2.071053

fit[[2]][1]
# [1] 5.006

fit[[2]][1,]
# Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#        5.006        3.428        1.462        0.246

fit$centers
#   Sepal.Length Sepal.Width Petal.Length Petal.Width
# 1     5.006000    3.428000     1.462000    0.246000
# 2     5.901613    2.748387     4.393548    1.433871
# 3     6.850000    3.073684     5.742105    2.071053

# code listing 20.1. Plotting the centroids from a k-means cluster analysis
set.seed(1234)
fit <- kmeans(iris[1:4], 3)
means <- fit$centers
library(reshape2)
dfm <- melt(means)
names(dfm) <- c("Cluster", "Measurement", "Centimeters")
dfm$Cluster <- factor(dfm$Cluster)
head(dfm)
#   Cluster  Measurement Centimeters
# 1       1 Sepal.Length    5.006000
# 2       2 Sepal.Length    5.901613
# 3       3 Sepal.Length    6.850000
# 4       1  Sepal.Width    3.428000
# 5       2  Sepal.Width    2.748387
# 6       3  Sepal.Width    3.073684

library(ggplot2)
ggplot(data = dfm,
       aes(x=Measurement, y=Centimeters, group=Cluster)) +
  geom_point(size=3, aes(shape=Cluster, color=Cluster)) +
  geom_line(size=1, aes(color=Cluster)) +
  ggtitle("Profiles for Iris Clusters")


# 20.1.2. Control structures

## for()

1:5
# [1] 1 2 3 4 5

for (i in 1:5) print(1:i)
# [1] 1
# [1] 1 2
# [1] 1 2 3
# [1] 1 2 3 4
# [1] 1 2 3 4 5

5:1
# [1] 5 4 3 2 1

for (i in 5:1) {
  print(1:i)
}
# [1] 1 2 3 4 5
# [1] 1 2 3 4
# [1] 1 2 3
# [1] 1 2
# [1] 1


## if() and else

## ifelse()
pvalues <- c(.0867, .0018, .0054, .1572, .0183, .5386)
results <- ifelse(pvalues < 0.05, "Significant", "Not Significant")
results
# [1] "Not Significant" "Significant"     "Significant"     "Not Significant"
# [5] "Significant"     "Not Significant"


# The same result can be accomplished with explicit loops using
pvalues <- c(.0867, .0018, .0054, .1572, .0183, .5386)
results <- vector(mode = "character", length = length(pvalues))
for (i in 1:length(pvalues)) {
  if (pvalues[i] < 0.05) results[i] <- "Significant"
  else results[i] <- "Not Significant"
}

results
# [1] "Not Significant" "Significant"     "Significant"     "Not Significant"
# [5] "Significant"     "Not Significant"

results[2]
# [1] "Significant"

# 20.1.3. Creating functions

## Function syntax

f <- function(x ,y, z=1) {
  result <- x + (2*y) + (3*z)
  return(result)
}

f(2, 3, 4)
# [1] 20

f(2,3)  # equal f(2, 3, 1)
# [1] 11

f(x=2, y=3) # equal f(x=2, y=3, z=1)
# [1] 11

f(z=4, y=2, 3)  # equal f(x=3, y=2, z=4)
# [1] 19
f(z=4, y=2, x=3)

args(f)
# function (x, y, z = 1) 
# NULL


## Object scope

x <- 2
y <- 3
z <- 4
f <- function(w) {
  z <- 2
  x <- w*y*z
  return(x)
}
f(x)
# [1] 12

x
# [1] 2

y
# [1] 3

z
# [1] 4


#============================================
# 20.2. Working with environments

# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

x <- 5
myenv <- new.env()
assign("x", "Homer", env=myenv)
ls()
# [1] "myenv" "x"
ls(myenv)
# [1] "x"
x
# [1] 5
get("x", env=myenv)
# [1] "Homer"

myenv$x
# [1] "Homer"

parent.env(myenv)
# <environment: R_GlobalEnv>
help("environment")

trim <- function(p){
  trimit <- function(x){
    n <- length(x)
    lo <- floor(n*p) + 1
    hi <- n + 1 - lo
    x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  trimit
}

x <- 1:10
x
# [1]  1  2  3  4  5  6  7  8  9 10

trim10pct <- trim(.1)
y <- trim10pct(x)
y
# [1] 2 3 4 5 6 7 8 9

trim20pct <- trim(.2)
y <- trim20pct(x)
y
# [1] 3 4 5 6 7 8

ls(environment(trim10pct))
# [1] "p"      "trimit"

get("p", env=environment(trim10pct))
# [1] 0.1

get("p", env=environment(trim20pct))
# [1] 0.2


makeFunction <- function(k) {
  f <- function(x) {
    print(x + k)
  }
}

g <- makeFunction(10)
g(4)
# [1] 14

k <- 2
g(5)
# [1] 15

ls(environment(g))
# [1] "f" "k"

environment(g)$k
# [1] 10


# 20.3.1 Generic functions
summary(women)
fit <- lm(weight ~ height, data = women)
summary(fit)

class(women)
class(fit)
# [1] "lm"

methods(summary)
# get code
summary.data.frame

# get code labeled with '*'
getAnywhere(summary.ggplot)

# code lisitng 20.2. An example of an arbitrary generic function
mymethod <- function(x, ...) UseMethod("mymethod")
mymethod.a <- function(x) print("Using A")
mymethod.b <- function(x) print("Using B")
mymethod.default <- function(x) print("Using Default")

x <- 1:5
y <- 6:10
z <- 10:15
class(x) <- "a"
class(y) <- "b"

mymethod(x)
# [1] "Using A"

mymethod(y)
# "Using B"

mymethod(z)
# "Using Default"

class(z) <- c("a", "b")
mymethod(z)
# [1] "Using A"

class(z) <- c("c", "a", "b")
mymethod(z)
# [1] "Using A"

#=======================================
# 20.3.2. Limitations of the S3 model
class(women) <- "lm"
summary(women)
# Error in if (p == 0) { : argument is of length zero

# 20.4. Writing efficient code

set.seed(1234)
mymatrix <- matrix(rnorm(10000000), ncol = 10)

accum <- function(x) {
  sums <- numeric(ncol(x))
  for (i in 1:ncol(x)) {
    for (j in 1:nrow(x)) {
      sums[i] <- sums[i] + x[j,i]
    }
  }
}

system.time(accum(mymatrix))
#  user  system elapsed 
# 0.732   0.002   0.736 

system.time(colSums(mymatrix))
#  user  system elapsed 
# 0.011   0.000   0.011


0.736/0.011

set.seed(1234)
k <- 100000
x <- rnorm(k)

y <- 0
system.time(for (i in 1:length(x)) y[i] <- x[i]^2)
#  user  system elapsed 
# 0.027   0.004   0.031 

y <- numeric(length=k)
system.time(for (i in 1:k) y[i] <- x[i]^2)
#  user  system elapsed 
# 0.008   0.001   0.008 


## Parallelization
# code listing 20.3. Parallelization with foreach and doParallel
library(foreach)
# install.packages("doParallel")
library(doParallel)
registerDoParallel(cores = 2) # Total Number of Cores: 2 from about this mac

eig <- function(n, p){
  x <- matrix(rnorm(100000), ncol = 100)
  r <- cor(x)
  eigen(r)$values
}
n <- 1000000
p <- 100
k <- 500

system.time(
  x <- foreach(i=1:k, .combine=rbind) %do% eig(n, p)
)
#  user  system elapsed 
# 7.612   0.416   8.064

system.time(
  x <- foreach(i=1:k, .combine = rbind) %dopar% eig(n, p)
)
#  user  system elapsed 
# 7.867   0.559   4.260 


#======================================================
# 20.5. Debugging
# 20.5.1. Common sources of errors
mtcars$Transmission <- factor(mtcars$a,
                              levels = c(1, 2),
                              labels = c("Automatic", "Manual"))
aov(mpg ~ Transmission, data = mtcars)
# Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
#   contrasts can be applied only to factors with 2 or more levels

head(mtcars[c("mpg", "Transmission")])
#                    mpg Transmission
# Mazda RX4         21.0    Automatic
# Mazda RX4 Wag     21.0    Automatic
# Datsun 710        22.8    Automatic
# Hornet 4 Drive    21.4         <NA>
# Hornet Sportabout 18.7         <NA>
# Valiant           18.1         <NA>

table(mtcars$Transmission)
# Automatic    Manual 
#        13         0 

# 20.5.2. Debugging tools
# code listing 20.4. A sample debugging session
args(mad)
# function (x, center = median(x), constant = 1.4826, na.rm = FALSE, 
#     low = FALSE, high = FALSE) 
# NULL

debug(mad)
mad(1:10)


# 20.5.3. Session options that support debugging
# code listing 20.5. Sample debugging session with recover()

f <- function(x, y) {
  z <- x + y
  g(z)
}

g <- function(x) {
  z <- round(x)
  h(z)
}

h <- function(x) {
  set.seed(1234)
  z <- rnorm(x)
  print(z)
}

options(error=recover)
f(2, 3)
# [1] -1.2070657  0.2774292  1.0844412 -2.3456977  0.4291247

f(2, -3)
# Error in rnorm(x) : invalid arguments
# 
# Enter a frame number, or 0 to exit   
# 
# 1: f(2, -3)
# 2: #3: g(z)
# 3: #3: h(z)
# 4: #3: rnorm(x)

# Selection: 4
# Called from: rnorm(x)
# Browse[1]> ls()
# [1] "mean" "n"    "sd"  
# Browse[1]> mean
# [1] 0
# Browse[1]> print(n)
# [1] -1
# Browse[1]> c
# Enter a frame number, or 0 to exit   
# 
# 1: f(2, -3)
# 2: #3: g(z)
# 3: #3: h(z)
# 4: #3: rnorm(x)
#   
# Selection: 3
# Called from: h(z)
# Browse[1]> ls()
# [1] "x"
# Browse[1]> x
# [1] -1
# Browse[1]> c
# 
# Enter a frame number, or 0 to exit   
# 
# 1: f(2, -3)
# 2: #3: g(z)
# 3: #3: h(z)
# 4: #3: rnorm(x)
#   
# Selection: 2
# Called from: g(z)
# Browse[1]> ls()
# [1] "x" "z"
# Browse[1]> x
# [1] -1
# Browse[1]> z
# [1] -1
# Browse[1]> c
# 
# Enter a frame number, or 0 to exit   
# 
# 1: f(2, -3)
# 2: #3: g(z)
# 3: #3: h(z)
# 4: #3: rnorm(x)
#   
# Selection: 1
# Called from: f(2, -3)
# Browse[1]> ls()
# [1] "x" "y" "z"
# Browse[1]> x
# [1] 2
# Browse[1]> y
# [1] -3
# Browse[1]> z
# [1] -1
# Browse[1]> c
# 
# Enter a frame number, or 0 to exit   
# 
# 1: f(2, -3)
# 2: #3: g(z)
# 3: #3: h(z)
# 4: #3: rnorm(x)
#   
# Selection: 0


options(error=NULL)
```
