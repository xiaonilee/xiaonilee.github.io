---
title: "Chapter 11. Intermediate graphs"
date: 2020-11-03
lastmod: 2020-11-03
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

Notebook of Reading Books: R in Action_Chapter 11.

<!--more-->

## This chapter covers

- Visualizing bivariate and multivariate relationships

- Working with scatter and line plots

- Understanding correlograms

- Using mosaic and association plots

### 11.1. Scatter plots

- with abline() and lowess()

- `id.method` option indicates that points will be identified interactively by mouse clicks, until the user selects Stop (via the Graphics or context-sensitive menu) or the Esc key. 

- Figure 11.1. Scatter plot of car mileage versus weight, with superimposed linear and lowess fit lines.
  
  ![fig111](fig111.png)

- Figure 11.2. Scatter plot with subgroups and separately estimated fit lines. 

  ![fig112](fig112.png)

#### 11.1.1. Scatter plot matrices

- Figure 11.3. Scatter plot matrix created by the `pairs()` function.

  ![fig113](fig113.png)

- Figure 11.3-1. Scatter plot matrix created by the `pairs()` function with **lower.panel = NULL**.

  ![fig1131](fig1131.png)

- Figure 11.3-2. Scatter plot matrix created by the `pairs()` function with **upper.panel = NULL**.
  
  ![fig1132](fig1132.png)

- Figure 11.4. Scatter plot matrix created with the `scatterplotMatrix()` function. 
  
  ![fig114](fig114.png)

- Figure 11.4-2. Scatter plot matrix produced with the `cpairs()` function in the `gclus` package. 
  - Variables closer to the principal diagonal are more highly correlated.
  
  ![fig1142](fig1142.png)

#### 11.1.2. High-density scatter plots

- Figure 11.5. Scatter plot with 10,000 observations and significant overlap of data points. 
  - Note that the overlap of data points makes it difficult to discern where the concentration of data is greatest.

  ![fig115](fig115.png)

- Figure 11.6. Scatterplot using `smoothScatter()` to plot **smoothed density estimates**. 
  - Densities are easy to read from the graph.

  ![fig116](fig116.png)

- Figure 11.7. Scatter plot using **hexagonal binning** to display the number of observations at each point. 
  - Data concentrations are easy to see and counts can be read from the legend.

  ![fig117](fig117.png)

- Figure 11.7-1. Scatter plot of 10,000 observations, where density is indicated by **color**. 
  - The data concentrations are easily discernable.

  ![fig1171](fig1171.png)

#### 11.1.3. 3D scatter plots

- scatterplot3d(x, y, z)

- Figure 11.8. `3D scatter plot` of miles per gallon, auto weight, and displacement.
  
  ![fig118](fig118.png)

- Figure 11.9. 3D scatter plot with vertical lines and shading.

  ![fig119](fig119.png)

- Figure 11.10. 3D scatter plot with vertical lines, shading, and **overlaid regression plane**.

  ![fig1110](fig1110.png)

#### 11.1.4 Spinning 3D Scatter Plots

- plot3d(x, y, z)

- Figure 11.11(up) and figure 11.11-1(down). Rotating 3D scatter plot produced by the `plot3d()` function in the `rgl` package.

  ![fig1111](fig1111.png)
  ![fig11111](fig11111.png)

- Figure 11.12 and figure 11.12-2. Spinning 3D scatter plot produced by the `scatter3d()` function in the `Rcmdr` package.

  ![fig1112](fig1112.png)
  ![fig11121](fig11121.png)

#### 11.1.5. Bubble plots

- Figure 11.13. Bubble plot of car weight versus mpg where point size is proportional to engine displacement.

  ![fig1113](fig1113.png)

### 11.2. Line charts

- Figure 11.14. **Comparison** of a scatter plot and a line plot.

  ![fig1114](fig1114.png)

- Table 11.1. Line chart options.

  ![tab111](tab111.png)

- Figure 11.15. type=options in the `plot()` and `lines()` functions.

  ![fig1115](fig1115.png)

- Figure 11.16. Line chart displaying the growth of five orange trees.

  ![fig1116](fig1116.png)

### 11.3. Correlograms

- with corrgram()

- Figure 11.17. Correlogram of the correlations among the variables in the `mtcars` data frame. 
  - Rows and columns have been reordered using **principal components analysis**.

  ![fig1117](fig1117.png)

- Figure 11.18. Correlogram of the correlations among the variables in the mtcars data frame. 
  - The lower triangle contains smoothed best fit lines and confidence ellipses; 
  - And the upper triangle contains scatter plots. 
  - The diagonal panel contains minimum and maximum values. 
  - Rows and columns have been reordered using principal components analysis.

  ![fig1118](fig1118.png)

- Figure 11.19. Correlogram of the correlations among the variables in the mtcars data frame. 
  - The lower triangle is shaded to represent the magnitude and direction of the correlations. 
  - The variables are plotted in their **original order**.

- figure 11.20. A Corrgram (or Horse) of a Different Color.

  ![fig1120](fig1120.png)

### 11.4. Mosaic plots

- Figure 11.21. `Mosaic` plot describing Titanic survivors by class, sex, and age.

  ![fig1121](fig1121.png)

Attach is the [Script](chapter11.R) of chapter11.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 11.1. Scatter plots
# code listing 11.1. A scatter plot with best fit lines
attach(mtcars)
plot(wt, mpg,
     main = "Basic Scatter plot of MPG vs. Weight",
     xlab="Car Weight (lbs/1000)",
     ylab="Miles Per Gallon ", pch=19)

abline(lm(mpg ~ wt), col="red", lwd=2, lty=1)
lines(lowess(wt, mpg), col="blue", lwd=2, lty=2) #figure 11.1

# figure 11.2
library(car)
scatterplot(mpg ~ wt | cyl, data=mtcars, lwd=2, span=0.75,
            main="Scatter Plot of MPG vs. Weight by # Cylinders",
            xlab="Weight of Car (lbs/1000)",
            ylab="Miles Per Gallon",
            legend.plot=TRUE,
            id=list(method="identify"),
            #labels=row.names(mtcars),
            boxplots="xy"
)


# 11.1.1. Scatter plot matrices
# figure 11.3
pairs(~mpg+disp+drat+wt, data = mtcars,
      main="Basic Scatter Plot Matrix")

# figure 11.3-1
pairs(~mpg+disp+drat+wt, data = mtcars,
      lower.panel = NULL,
      main="Basic Scatter Plot Matrix")

# figure 11.3-2
pairs(~mpg+disp+drat+wt, data = mtcars,
      upper.panel = NULL,
      main="Basic Scatter Plot Matrix")

# figure 11.4
library(car)
scatterplotMatrix(~ mpg + disp + drat + wt, data = mtcars, spead=F,
                  smoother.args=list(lty=2),
                  main="Scatter Plot Matrix via car Package")


cor(mtcars[c("mpg", "wt", "disp", "drat")])
#             mpg         wt       disp       drat
# mpg   1.0000000 -0.8676594 -0.8475514  0.6811719
# wt   -0.8676594  1.0000000  0.8879799 -0.7124406
# disp -0.8475514  0.8879799  1.0000000 -0.7102139
# drat  0.6811719 -0.7124406 -0.7102139  1.0000000
# There were 50 or more warnings (use warnings() to see the first 50)

# code listing 11.2. Scatter plot matrix produced with the gclus package
# install.packages("gclus")
library(gclus)
colnames(mtcars)
# [1] "mpg"  "cyl"  "disp" "hp"   "drat" "wt"   "qsec" "vs"   "am"   "gear" "carb"

str(mtcars)
mydata <- mtcars[c(1, 3, 5, 6)]
mydata.corr <- abs(cor(mydata))

mycolors <- dmat.color(mydata.corr)

myorder <- order.single(mydata.corr)

# figure 11.4-2
cpairs(mydata,
       myorder,
       panel.colors=mycolors,
       gap=.5,
       main="Variables Ordered and Colored by Correlation"
)

# 11.1.2. High-density scatter plots

set.seed(1234)

n <- 10000
c1 <- matrix(rnorm(n, mean=0, sd=.5), ncol=2)
c2 <- matrix(rnorm(n, mean=3, sd=2), ncol=2)
mydata <- rbind(c1, c2)
mydata <- as.data.frame(mydata)
names(mydata) <- c("x", "y")

# figure 11.5
with(mydata,
     plot(x, y, pch=19, main="Scatter Plot with 10,000 Observations"))

# figure 11.6
with(mydata,
     smoothScatter(x, y, main="Scatterplot Colored by Smoothed Densities"))

# figure 11.7
# install.packages("hexbin")
library(hexbin)
with(mydata, {
  bin <- hexbin(x, y, xbins=50)
  plot(bin, main="Hexagonal Binning with 10,000 Observations")
})


# figure 11.7-1
# install.packages("IDPmisc")
library(IDPmisc)
with(mydata,
     iplot(x, y, main="Image Scatter Plot with Color Indicating Density"))


# 11.1.3. 3D scatter plots
# figure 11.8
# install.packages("scatterplot3d")
library(scatterplot3d)
attach(mtcars)

scatterplot3d(wt, disp, mpg,
              main="Basic 3D Scatter Plot")


# figure 11.9
library(scatterplot3d)
attach(mtcars)
scatterplot3d(wt, disp, mpg,
              pch=16,
              highlight.3d = T,
              type = "h",
              main="3D Scatter Plot with Vertical Lines")

# figure 11.10
library(scatterplot3d)
attach(mtcars)
s3d <- scatterplot3d(wt, disp, mpg,
                     pch=16,
                     highlight.3d = T,
                     type = "h",
                     main="3D Scatter Plot with Vertical Lines")
fit <- lm(mpg ~ wt+disp)
s3d$plane3d(fit)

# figure 11.11 and 11.11-1
# continue figure 11.10
# install.packages("rgl")
library(rgl)
attach(mtcars)
plot3d(wt, disp, mpg, col = "red", size=5)

# figure 11.12 and figure 11.12-1
# install.packages("Rcmdr")
library(Rcmdr)
attach(mtcars)
scatter3d(wt, disp, mpg)


# 11.1.5. Bubble plots
attach(mtcars)
r <- sqrt(disp/pi)
symbols(wt, mpg, circles = r, inches = 0.30,
        fg="white", bg="lightblue",
        main="Bubble Plot with point size proportional to displacement",
        ylab = "Miles Per Gallon",
        xlab = "Weight of Car (lbs/1000)")
text(wt, mpg, rownames(mtcars), cex=0.6)
detach(mtcars)

# 11.2. Line charts
# code listing 11.3. Creating side-by-side scatter and line plots
opar <- par(no.readonly = T)
par(mfrow=c(1,2))
table(Orange$Tree)
t1 <- subset(Orange, Tree==1)
plot(t1$age, t1$circumference,
     xlab = "Age (days)",
     ylab = "Circumference (mm)",
     main = "Orange Tree 1 Growth")
plot(t1$age, t1$circumference,
     xlab = "Age (days)",
     ylab = "Circumference (mm)",
     main = "Orange Tree 1 Growth",
     type = "b")
par(opar)

# figure 11.15
opar <- par(no.readonly = T)
par(mfrow=c(2,4))

linetpye <- c("p", "l", "o", "b", "c", "s", "S", "h")

t1 <- subset(Orange, Tree==1)

for (i in linetpye){
  mainstr <- paste("type = ", '"', i, '"')
  plot(t1$age, t1$circumference,
       xlab = "Age (days)",
       ylab = "Circumference (mm)",
       type = "n",
       main = mainstr)
  lines(t1$age, t1$circumference, type = i)
}

par(opar)


# code listing 11.4. Line chart displaying the growth of five orange trees over time
# figure 11.16
class(Orange$Tree)
Orange$Tree <- as.numeric(Orange$Tree)
class(Orange$Tree)

ntrees <- max(Orange$Tree)

xrange <- range(Orange$age)
yrange <- range(Orange$circumference)

plot(xrange, yrange,
     type = "n",
     xlab = "Age (days)",
     ylab = "Circumference (mm)"
     )

colors <- rainbow(ntrees)
linetype <- c(1:ntrees)
plotchar <- seq(18, 18+ntrees, 1)

for (i in 1:ntrees){
    tree <- subset(Orange, Tree==i)
    lines(tree$age, tree$circumference,
          type="b",
          lwd=2,
          lty=linetype[i],
          col=colors[i],
          pch=plotchar[i]
    )
}

title("Tree Growth", "example of line plot")

legend(xrange[1], yrange[2],
       1:ntrees,
       cex = 0.8,
       col = colors,
       pch = plotchar,
       lty=linetype,
       title = "Tree"
)

# 11.3. Correlograms
options(digits = 2)
cor(mtcars)
#        mpg   cyl  disp    hp   drat    wt   qsec    vs     am  gear   carb
# mpg   1.00 -0.85 -0.85 -0.78  0.681 -0.87  0.419  0.66  0.600  0.48 -0.551
# cyl  -0.85  1.00  0.90  0.83 -0.700  0.78 -0.591 -0.81 -0.523 -0.49  0.527
# disp -0.85  0.90  1.00  0.79 -0.710  0.89 -0.434 -0.71 -0.591 -0.56  0.395
# hp   -0.78  0.83  0.79  1.00 -0.449  0.66 -0.708 -0.72 -0.243 -0.13  0.750
# drat  0.68 -0.70 -0.71 -0.45  1.000 -0.71  0.091  0.44  0.713  0.70 -0.091
# wt   -0.87  0.78  0.89  0.66 -0.712  1.00 -0.175 -0.55 -0.692 -0.58  0.428
# qsec  0.42 -0.59 -0.43 -0.71  0.091 -0.17  1.000  0.74 -0.230 -0.21 -0.656
# vs    0.66 -0.81 -0.71 -0.72  0.440 -0.55  0.745  1.00  0.168  0.21 -0.570
# am    0.60 -0.52 -0.59 -0.24  0.713 -0.69 -0.230  0.17  1.000  0.79  0.058
# gear  0.48 -0.49 -0.56 -0.13  0.700 -0.58 -0.213  0.21  0.794  1.00  0.274
# carb -0.55  0.53  0.39  0.75 -0.091  0.43 -0.656 -0.57  0.058  0.27  1.000

# install.packages("corrgram")
library(corrgram)
corrgram(mtcars, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Correlogram of mtcars intercorrelations")  # figure 11.17



# a second example
# figure 11.18
library(corrgram)
corrgram(mtcars, order=TRUE, lower.panel=panel.ellipse,
         upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.minmax,
         main="Correlogram of mtcars data using scatter plots and ellipses")

# figure 11.19 with originalß order
library(corrgram)
corrgram(mtcars, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="Car Mileage Data (unsorted)")

# figure 11.20
library(corrgram)
col.corrgram <- function(ncol){
  colorRampPalette(c("darkgoldenrod4", "burlywood1",
                     "darkkhaki", "darkgreen"))(ncol)}
corrgram(mtcars, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="A Corrgram (or Horse) of a Different Color")



# 11.4. Mosaic plots
ftable(Titanic)
#                    Survived  No Yes
# Class Sex    Age                   
# 1st   Male   Child            0   5
#              Adult          118  57
#       Female Child            0   1
#              Adult            4 140
# 2nd   Male   Child            0  11
#              Adult          154  14
#       Female Child            0  13
#              Adult           13  80
# 3rd   Male   Child           35  13
#              Adult          387  75
#       Female Child           17  14
#              Adult           89  76
# Crew  Male   Child            0   0
#              Adult          670 192
#       Female Child            0   0
#              Adult            3  20

# figure 11.21
library(vcd)
mosaic(Titanic, shade=TRUE, legend=TRUE)

# altinative code
# mosaic(~Class+Sex+Age+Survived, data=Titanic, shade=TRUE, legend=TRUE)
```
