---
title: "Chapter 5. Advanced data management"
date: 2020-10-23
lastmod: 2020-10-28
draft: false
tags: ["R", "R in Action", "Bioinformatics", "Book"]
categories: ["R", "Bioinformatics"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

Notebook of Reading Books: R in Action_Chapter 5.

<!--more-->

### This chapter covers

- Mathematical and statistical functions
  
- Character functions
  
- Looping and conditional execution
  
- User-written functions
  
- Ways to aggregate and reshape data

### Notes about some useful functions

- 5.2.2. Statistical functions

  - scale()

  - tansform()

- 5.2.3. Probability functions

  - runif()

  - set.sead()

  - mvrnorm()

- 5.2.4. Character functions

  - nchar()

  - substr()

  - grep()

  - sub()

  - strsplit()

  - paste()

  - toupper() and tolower()

  - seq(from, to, by)

  - rep()

  - cut()

  - pretty()

  - cat()

- 5.2.6. Applying functions to matrices and data frames

  - apply(x, margin, FUN, ...)
  
- 5.4. Control flow

  - for (var in seq) statement

  - if (cond) statement [else statement2]

  - ifelse(cond, statement1, statement2)

  - switch(expr, ...)

- 5.6.2. Aggregating data

  - aggregate(x, by, FUN)

- 5.6.3. The reshape package, versatility

  - melt()

  - cast() and dcast()

Attach is the [Script](chapter5.R) of chapter5.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 5.1. A data management challenge

# 5.2. Numerical and character functions

# 5.2.1. Mathematical functions

# 5.2.2. Statistical functions

newdata <- scale(mydata)*SD + M

newdata <- transform(mydata, myvar = scale(myvar)*10+50)

# 5.2.3. Probability functions

# Setting the Seed for Random Number Generation
# Listing 5.2. Generating pseudo-random numbers from a uniform distribution
runif(5)
runif(5)
runif(5)
runif(5)

# Generating Multivariate Normal Data
# code listing 5.3
library(MASS)
options(digits = 3)
set.seed(1234)

mean <- c(230.7, 146.7, 3.6)
sigma <- matrix(c(15360.8, 6721.2, -47.1,
                  6721.2, 4700.9, -16.5,
                  -47.1, -16.5, 0.3), nrow = 3, ncol = 3)
sigma

mydata <- mvrnorm(500, mean, sigma)
mydata <- as.data.frame(mydata)
names(mydata) <- c("y", "x1", "x2")
head(mydata)
dim(mydata)

# 5.2.4. Character functions
x <- c("ab", "cde", "fghij")
nchar(x)
# [1] 2 3 5

substr(x, 2, 4)
# [1] "b"   "de"  "ghi"

substr(x, 2, 4) <- "2222"
# [1] "a2"    "c22"   "f222j"

grep(pattern, x, ignore.case=FALSE, fixed=FALSE)
grep("A", c("b","A","c"), fixed=TRUE)
# [1] 2

sub(pattern, replacement, x, ignore.case=FALSE, fixed=FALSE)
sub("\\s",".","Hello There")
# [1] "Hello.There"

strsplit(x, split, fixed=FALSE)
y <- strsplit("abc", "")
y
# [[1]]
# [1] "a" "b" "c"

paste(..., sep="")
paste("x", 1:3, sep="")
# c("xl", "x2", "x3")
paste("x",1:3,sep="M")
# c("xMl","xM2" "xM3")

toupper(x)
tolower(x)

# 5.2.5. Other useful functions
length(x)
length(c(2,5,6,9))

seq(from,to,by)
seq(1,10,2)

rep(x,n)
rep(1:3,2)
# [1] 1 2 3 1 2 3

cut(x,n)

pretty(x,n)

cat(..., file = "myfile", append = FALSE)
firstname <- c("Jane")
cat("Hello",firstname, "\n") # equal cat("HEllo", "Jane", "\n")

?Quotes

name <- "Bob"
cat( "Hello", name, "\b.\n", "Isn\'t R", "\t", "GREAT?\n")
cat( "Hello", name, ".\n", "Isn\'t R", "\t", "GREAT?\n")

# 5.2.6. Applying functions to matrices and data frames

apply(x, margin, FUN, ...)

# code listing 5.5
mydata <- matrix(rnorm(30), nrow = 6)
mydata

apply(mydata, 1, mean)
apply(mydata, 2, mean)

# Calculate trimmer column means
apply(mydata, 2, mean, trim=0.2) # trim topmost 20% and lowest20%


# 5.3. A solution for our data management challenge
# code listing 5.6
options(digits = 2)
Student <- c("John Davis", "Angela Williams", "Bullwinkle Moose",
             "David Jones", "Janice Markhammer", "Cheryl Cushing",
             "Reuven Ytzrhak", "Greg Knox", "Joel England",
             "Mary Rayburn")
Math <- c(502, 600, 412, 358, 495, 512, 410, 625, 573, 522)
Science <- c(95, 99, 80, 82, 75, 85, 80, 95, 89, 86)
English <- c(25, 22, 18, 15, 20, 28, 15, 30, 27, 18)

roster <- data.frame(Student, Math, Science, English,
                     stringsAsFactors = F)

head(roster)
# Obtain performance scores
z <- scale(roster[,2:4])
z
score <- apply(z, 1, mean)
roster <- cbind(roster, score)

# Grade students
y <- quantile(score, c(0.8, 0.6, 0.4, 0.2))

roster$grade[score >= y[1]] <- "A"
roster$grade[score < y[1] & score >= y[2]] <- "B"
roster$grade[score < y[2] & score >= y[3]] <- "C"
roster$grade[score < y[3] & score >= y[4]] <- "D"
roster$grade[score < y[4]] <- "F"

# Extract last and first names
name <- strsplit((roster$Student), " ")
Lastname <- sapply(name, "[", 2)
Firstname <- sapply(name, "[", 1)
?"[" # Extract or Replace Parts of an Object


roster <- cbind(Firstname, Lastname, roster[, -1])
head(roster)

# Sort by last and first names
roster <- roster[order(Lastname,Firstname),]
roster

# 5.4. Control flow
# 5.4.1. Repetition and looping

# for (var in seq) statement

for (i in 1:10) print("Hello")

# while (cond) statement
i <- 10
while (i > 0) {print("Hello"); i <- i - 1}

# 5.4.2. Conditional execution
# if (cond) statement
# if (cond) statement else statement2

if (is.character(grade)) grade <- as.factor(grade)
if (!is.character(grade)) grade <- as.factor(grade) else print("Grade already
                                                              is a factor")

# ifelse(cond, statement1, statement2)
ifelse(score > 0.5, print("Passed"), print("Failed"))
outcome <- ifelse(score > 0.5, "Passed", "Failed")

# switch(expr, ...)
# code listing 5.7
feelings <- c("sad", "afraid")
for (i in feelings)
  print(
    switch(i,
           happy  = "I am glad you are happy",
           afraid = "There is nothing to fear",
           sad    = "Cheer up",
           angry  = "Calm down now"
    )
  )

# 5.5. User-written functions

# The structure of a function looks like this:
myfunction <- function(arg1,arg2,...) {
  statements
  return(object)
}

# code listing 5.8
mystats <- function(x, parametric=TRUE, print=FALSE) {
  if (parametric) {
    center <- mean(x); spread <- sd(x)
  } else {
    center <- median(x); spread <- mad(x)
  }
  if (print & parametric) {
    cat("Mean=", center, "\n", "SD=", spread, "\n")
  } else if (print & !parametric) {
    cat("Median=", center, "\n", "MAD=", spread, "\n")
  }
  result <- list(center=center, spread=spread)
  return(result)
}


set.seed(1234)
x <- rnorm(500)

y <- mystats(x) # if (parametric)
y <- mystats(x, parametric = F, print = T) # else if (print & !parametric)
y <- mystats(x, parametric = T, print = T) # if (print & parametric)

# a user-written function that uses the switch construct
mydate <- function(type="long") {
  switch(type,
         long =  format(Sys.time(), "%A %B %d %Y"),
         short = format(Sys.time(), "%m-%d-%y"),
         cat(type, "is not a recognized type\n")
  )
}
mydate("long")
mydate("short")
mydate("middle")

# 5.6. Aggregation and restructuring
# 5.6.1. Transpose

# code listing 5.9
cars <- mtcars[1:5,1:4]
head(cars)
t(cars)

# 5.6.2. Aggregating data
# aggregate(x, by, FUN)

# code listing 5.10
options(digits=3)
attach(mtcars)
aggdata <-aggregate(mtcars, by=list(Group.cyl=cyl,Group.gear=gear), FUN=mean, na.rm=TRUE)
head(aggdata)

# 5.6.3. The reshape package, versatility/万能特性

ID <- c(1, 1, 2, 2)
Time <- c(1, 2, 1, 2)
X1 <- c(5, 3, 6, 2)
X2 <- c(6, 5, 1, 4)
mydata <- data.frame(ID, Time, X1, X2)
head(mydata)

# Melting
library(reshape2)
md <- melt(mydata, id=(c("ID", "Time")))
md

# Casting
newdata <- cast(md, formula, FUN)

dcast(md, ID+Time~variable)
dcast(md, ID+variable~Time)
dcast(md, ID~variable+Time)

dcast(md, ID~variable,mean)
dcast(md,ID~Time,mean)


# 5.7 Summary
```
