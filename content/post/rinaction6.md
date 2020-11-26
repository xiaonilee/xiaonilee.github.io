---
title: "Chapter 6. Basic graphs"
date: 2020-10-24
lastmod: 2020-10-28
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

Notebook of Reading Books: R in Action_Chapter 6.

<!--more-->

### This chapter covers

- Bar, box, and dot plots
  
- Pie and fan charts
  
- Histograms and kernel density plots

### Note about pie charts

- **Figure for Pie charts with code Listing 6.5.**

![piecharts](piecharts.png)

### Note about histogram charts

- **Figure for histogram charts with code Listing 6.6.**

- The surrounding box is produced by the ***box()*** function.

![hist6.6](hist6.6.png)

### Note about density plots

- **Figure for density plots with code Listing 6.7 and 6.8.**

- The ***locator(1)*** option indicates that you’ll place the legend interactively by clicking on the graph where you want the legend to appear.

![density](densityplot.png)

### Note about box plots

- **Figure for box plots with code listing 6.9.**

![boxplot](boxplot1.png)

### Note about violin plots

- **Figure for violin plots with code listing 6.10.**

- A `violin plot` is a **combination** of a ***box plot*** and a ***kernel density plot***.

![violinplot](violinplot.png)

### Note about dot charts

- **Figure for dot charts with code listing 6.11.**

- Dot plot of mpg for car models **grouped**, **sorted**, and **colored**.

![dotchart](dotchart.png)

Attach is the [Script](chapter6.R) of chapter6.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)


# 6.1. Bar plots
library(vcd)
head(Arthritis)
dim(Arthritis)

# 6.1.1. Simple bar plots
counts <- table(Arthritis$Improved)
counts

# code listing 6.1
barplot(counts,
        main = "Simple Bar Plot",
        xlab = "Improvement",
        ylab = "Frequency")

barplot(counts,
        main = "Horizontal Bar Plot",
        ylab = "Improvement",
        xlab = "Frequency",
        horiz = T)

# 6.1.2. Stacked and grouped bar plots
counts <- table(Arthritis$Improved, Arthritis$Treatment)
counts

# code listing 6.2
barplot(counts,
        main = "Stacked Bar Plot",
        xlab = "Treatment",
        ylab = "Frequency",
        col = c("red", "yellow", "green"),
        legend=rownames(counts))

barplot(counts,
        main = "Grouped Bar Plot",
        xlab = "Treatment",
        ylab = "Frequency",
        col = c("red", "yellow", "green"),
        legend=rownames(counts), beside = T)

# 6.1.3. Mean bar plots
# code listing 6.3
states <- data.frame(state.region, state.x77)
head(states)

means <- aggregate(states$Illiteracy, by=list(state.region), FUN=mean)
means

means <- means[order(means$x),]

barplot(means$x, names.arg = means$Group.1, title("Mean Illiteracy Rate"))

# 6.1.4. Tweaking bar plots
# code listing 6.4
par(mar=c(5,8,4,2))
par(las=2)
counts <- table(Arthritis$Improved)

barplot(counts,
        main="Treatment Outcome",
        horiz=TRUE, cex.names=0.8,
        names.arg=c("No Improvement", "Some Improvement",
                    "Marked Improvement"))

# 6.1.5. Spinograms spine()

library(vcd)
attach(Arthritis)
counts <- table(Treatment, Improved)
spine(counts, main="Spinogram Example")
detach(Arthritis)

# 6.2. Pie charts

pie(x,labels)  # pie charts
fan.plot()     # fan plot, a variation of the pie chart

# code Listing 6.5. Pie charts
par(mfrow=c(2,2))
slices <- c(10, 12, 4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")

pie(slices, labels = lbls,
    main = "Simple Pie Chart")

pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls, " ", pct, "%", sep = "")
pie(slices, labels = lbls2, col = rainbow(length(lbls2)),
    main = "Pie Chart with Percentages")

library(plotrix)
pie3D(slices, labels=lbls, explode=0.1,
      main="3D Pie Chart")

mytable <- table(state.region)
lbls3 <- paste(names(mytable), "\n", mytable, sep = "")
pie(mytable, labels = lbls3,
    main = "Pie Chart from a Table\n (with sample sizes)")


# 6.3. Histograms

# code listing 6.6. Histograms

par(mfrow=c(2,2))

hist(mtcars$mpg)

hist(mtcars$mpg,
     breaks = 12,
     col = "red",
     xlab = "Miles Per Gallon",
     main = "Colored histogram with 12 bins")

hist(mtcars$mpg,
     freq = F,
     breaks = 12,
     col = "red",
     xlab = "Miles Per Gallon",
     main = "Histogram, rug plot, density curve")
# rug plot
rug(jitter(mtcars$mpg)) # a one-dimensional representation of the actual data values

lines(density(mtcars$mpg), col="blue", lwd=2)

x <- mtcars$mpg
h <- hist(x,
          breaks = 12,
          col = "red",
          xlab = "Miles Per Gallon",
          main = "Histogram with normal curve and box")
xfit <- seq(min(x), max(x), length=40)
yfit <- dnorm(xfit, mean = mean(x), sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)
box() # produce surrounding box

# 6.4. Kernel density plots

plot(density(x))

# if overlap  with an existing graph
lines()

# code listing 6.7. Kernel density plots
par(mfrow=c(2, 2))
d <- density(mtcars$mpg)
plot(d)

d <- density(mtcars$mpg)
plot(d, main = "Kernel Density of Miles Per Gallon")
polygon(d, col = "red", border = "blue")
rug(mtcars$mpg, col = "brown")

# code listing 6.8. Comparative kernel density plots
par(lwd=2)
library(sm)
attach(mtcars)

cyl.f <- factor(cyl, levels = c(4, 6, 8),
                labels = c("4 cylinder", "6 cylinder",
                           "8 cylinder"))

sm.density.compare(mpg, cyl, xlab="Miles Per Gallon")
title(main = "MPG Distribution by Car Cylinders")

colfill <- c(2:(1+length(levels(cyl.f))))

# locator(1), indicates that you’ll place the legend interactively
# by clicking on the graph where you want the legend to appear.
legend(locator(1), levels(cyl.f), fill = colfill)

detach(mtcars)

# 6.5. Box plots

boxplot.stats(mtcars$mpg)

## 6.5.1. Using parallel box plots to compare groups
table(mtcars$cyl)
table(mtcars$mpg)
summary(mtcars$mpg)
table(mtcars$mpg,  mtcars$cyl)
table(mtcars$mpg, mtcars$cyl)

# code listing 6.9. Box plots for two crossed factors

mtcars$cyl.f <- factor(mtcars$cyl,
                       levels = c(4, 6, 8),
                       labels = c("4", "6", "8"))

mtcars$am.f <- factor(mtcars$am,
                      levels = c(0, 1),
                      labels = c("auto", "standard"))

boxplot(mpg ~ am.f * cyl.f,
        data = mtcars,
        varwidth=T,
        col=c("gold", "darkgreen"),
        main="MPG Distribution by Auto Type",
        xlab = "Auto Type")

# 6.5.2. Violin plots

# A violin plot is a combination of a box plot and a kernel density plot.
vioplot(x1, x2, ..., names=, col=)

# code listing 6.10. Violin plots
library(vioplot)
x1 <- mtcars$mpg[mtcars$cyl == 4]
x2 <- mtcars$mpg[mtcars$cyl == 6]
x3 <- mtcars$mpg[mtcars$cyl == 8]
vioplot(x1, x2, x3,
        names=c("4 cyl", "6 cyl", "8 cyl"),
        col="gold")
title("Violin Plots of Miles Per Gallon")

# 6.6. Dot plots
# code listing 6.11.Dot plot grouped, sorted, and colored
x <- mtcars[order(mtcars$mpg),]
x$cyl <- factor(x$cyl)
x$color[x$cyl==4] <- "red"
x$color[x$cyl==6] <- "blue"
x$color[x$cyl==8] <- "darkgreen"
dotchart(x$mpg,
         labels = row.names(x),
         cex=.7,
         groups = x$cyl,
         gcolor = "black",
         color = x$color,
         pch=19,
         main = "Gas Mileage for Car Models\ngrouped by cylinder",
         xlab = "Miles Per Gallon")
```
