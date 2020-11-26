---
title: "Chapter 23. Advanced graphics with the lattice package"
date: 2020-11-23
lastmod: 2020-11-23
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

Notebook of Reading Books: R in Action_Chapter 23.

<!--more-->

## This chapter covers

- An introduction to the lattice package

- Grouping and conditioning

- Adding information with panel functions

- Customizing a lattice graph’s appearance

### 23.1. The lattice package

- Figure 23.1. Trellis graph of singer heights by voice part.

  ![fig231](fig231.png)

- Table 23.1. Graph types and corresponding functions in the ***lattice*** package.

  ![tab231](tab231.png)

### 23.2. Conditioning variables

- Typically, conditioning variables are factors. 

- If It's a continuous variable:
  - `cut()`
  - **myshingle**

```r
myshingle <- equal.count(x, number=n, overlap=proportion)
```

- Figure 23.2. Trellis plot of miles per gallon vs. car weight conditioned on engine displacement.
  - engine displacement is a continuous variable
    - it has been converted to three non-overlapping shingles with equal numbers of observations.

  ![fig232](fig232.png)

### 23.3. Panel functions

- Figure 23.3. Trellis plot of miles per gallon vs. car weight conditioned on engine displacement. 
  - A custom panel function has been used to add regression lines, rug plots, and grid lines.

  ![fig233](fig233.png)

- Figure 23.4. Trellis graph of miles per gallon vs. engine displacement conditioned on transmission type. 
  - Smoothed lines (loess), grids, and group mean levels have been added.

  ![fig234](fig234.png)

### 23.4. Grouping variables

- Figure 23.5. Kernel-density plots for miles per gallon grouped by transmission type. 
  - Jittered points are provided on the horizontal axis.

  ![fig235](fig235.png)

- Figure 23.6. Kernel-density plots for miles per gallon grouped by transmission type. 
  - Graphical parameters have been modified, and a customized legend has been added. 
  - The custom legend specifies color, shape, line type, character size, and title.

  ![fig236](fig236.png)

- Figure 23.7. xyplot showing the impact of ambient carbon dioxide concentrations on carbon dioxide uptake for 12 plants in two treatment conditions and two types. 
  - Plant is the group variable, and Treatment and Type are the conditioning variables.

  ![fig237](fig237.png)

### 23.5. Graphic parameters

- The default of `show.settings()`

  ![fig23.0](fig23.0.png)

- The changes of `show.settings()`
  
  ![fig23.0-1](fig23.0-1.png)

### 23.6. Customizing plot strips

- Figure 23.8. A trellis graph with a customized strip (light grey background, with a smaller, italicized font).

  ![fig238](fig238.png)

### 23.7. Page arrangement

- Figure 23.9. Using the `split` option to combine graphs.

  ![fig239](fig239.png)

- Figure 23.10. Using the `position` option to combine graphs with greater precision.

  ![fig2310](fig2310.png)

Attach is the [Script](chapter23.R) of chapter23.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 23.1. The lattice package
library(lattice)
histogram(~height | voice.part, data = singer,
          main = "Distribution of Heights by Voice Pitch",
          xlab = "Height (inches)")

# code listing 23.1. Lattice plot examples
library(lattice)
attach(mtcars)

gear <- factor(gear, levels = c(3, 4, 5),
               labels = c("3 gears", "4 gears", "5 gears"))

cyl <- factor(cyl, levels = c(4, 6, 8),
              labels = c("4 cylinders", "6 cylinders", "8 cylinders"))

densityplot(~mpg,
            main = "Density Plot",
            xlab = "Miles per Gallon")

densityplot(~mpg | cyl,
            main = "Density Plot by Number of Cylinders",
            xlab = "Miles per Gallon")

bwplot(cyl ~ mpg | gear,
       main = "Box Plots by Cylinders and Gears",
       xlab = "Miles per Gallon", ylab = "Cylinders")

xyplot(mpg ~ wt | cyl * gear,
       main = "Scatter Plots by Cylinders and Gears",
       xlab = "Car Weight", ylab = "Miles per Gallon")

cloud(mpg ~ wt * qsec | cyl,
      main = "3D Scatter Plots by Cylinders")

dotplot(cyl ~ mpg | gear,
        main = "Dot Plots by Number of Gears and Cylinders",
        xlab = "Miles per Gallon")

splom(mtcars[c(1, 3, 4, 5, 6)],
      main = "Scatter Plot Matrix for mtcars Data")

detach(mtcars)

# 23.2. Conditioning variables

# figure 23-2
displacement <- equal.count(mtcars$disp, number=3, overlap=0)
xyplot(mpg~wt|displacement, data=mtcars,
       main = "Miles per Gallon vs. Weight by Engine Displacement",
       xlab = "Weight", ylab = "Miles per Gallon",
       layout=c(3, 1), aspect=1.5)


# 23.3. Panel functions
# code listing 23.2. xyplot with custom panel function
# figure 23.3
library(lattice)
displacement <- equal.count(mtcars$disp, number=3, overlap=0)

mypanel <- function(x, y) {
  panel.xyplot(x, y, pch = 19)
  panel.rug(x, y)
  panel.grid(h=-1, v=-1)
  panel.lmline(x, y, col="red", lwd=1, lty=2)
}

xyplot(mpg~wt | displacement, data = mtcars,
       layout=c(3, 1),
       aspect = 1.5,
       main = "Miles per Gallon vs. Weight by Engine Displacement",
       xlab = "Weight",
       ylab = "Miles per Gallon",
       panel = mypanel)


# code listing 23.3. xyplot with a custom panel function and additional options
# figure 23.4
library(lattice)
mtcars$transmission <- factor(mtcars$am, levels = c(0, 1),
                              labels = c("Automatic", "Manual"))

panel.smoother <- function(x, y) {
  panel.grid(h=-1, v=-1)
  panel.xyplot(x, y)
  panel.loess(x, y)
  panel.abline(h=mean(y), lwd=2, lty=2, col = "darkgreen")
}

xyplot(mpg~disp | transmission, data = mtcars,
       scales = list(cex=.8, col="red"),
       panel = panel.smoother,
       xlab = "Displacement", ylab = "Miles per Gallon",
       main = "MPG vs Displacement by Transmission Type",
       sub = "Dotted lines are Group Means", aspect = 1)

#-----------------------------------------------------------
# Subsetting Data                                          ｜
## Selection using the Subset Function                     ｜ 
sub <- subset(mtcars, transmission == "Manual")           #｜
mean(sub$mpg)                                             #｜
#                                                          ｜
# Selecting Observations                                   ｜
mean(mtcars[which(mtcars$transmission=='Manual'),]$mpg)   #｜
# [1] 24.39231                                            #｜
mean(mtcars[which(mtcars$transmission=='Automatic'),]$mpg)#｜
# [1] 17.14737                                             ｜
#----------------------------------------------------------｜

# 23.4. Grouping variables
# figure 23.5
library(lattice)
mtcars$transmission <- factor(mtcars$am, levels = c(0, 1),
                              labels = c("Automatic", "Manual"))

densityplot(~mpg, data = mtcars,
            group = transmission,
            main = "MPG Distribution by Transmission Type",
            xlab = "Miles per Gallon",
            auto.key = TRUE)

# code listing 23.4. Kernel-density plot with a group variable and customized legend
# figure 23.6
library(lattice)
mtcars$transmission <- factor(mtcars$am, levels = c(0, 1),
                              labels = c("Automatic", "Manual"))

# Color, line and point specifications
colors <- c("red", "blue")
lines <- c(1, 2)
points <- c(16, 17)

# Legend customization
key.trans <- list(title="Transmission",
                  space="bottom", columns=2,
                  text=list(levels(mtcars$transmission)),
                  points=list(pch=points, col=colors),
                  lines=list(col=colors, lty=lines),
                  cex.title=1, cex=.9)
# Density plot
densityplot(~mpg, data = mtcars,
            group = transmission,
            main = "MPG Distribution by Transmission Type",
            xlab = "Miles per Gallon",
            pch=points, lty=lines, col=colors,
            lwd=2, jitter=.005,
            key=key.trans)

# code listing 23.5. xyplot with group and conditioning variables and customized legend
# figure 23.7

head(CO2)
library(lattice)
colors <- "darkgreen"
symbols <- c(1:12)
linetype <- c(1:3)

key.species <- list(title="Plant",
                    space="right",
                    text=list(levels(CO2$Plant)),
                    points=list(pch=symbols, col=colors))

xyplot(uptake ~ conc | Type*Treatment, data = CO2,
       group = Plant,
       type="o",
       pch=symbols, col=colors, lty=linetype,
       main="Carbon Dioxide Uptake\nin Grass Plants",
       ylab = expression(paste("Uptake ",
                               bgroup("(", italic(frac("umol", "m"^2)), ")"))),
       xlab = expression(paste("Concentration ",
                               bgroup("(", italic(frac(mL, L)), ")"))),
       sub = "Grass Species: Eninochloa crus-galli",
       key=key.species)

# 23.5. Graphic parameters

show.settings() # view the current defaults

mysettings <- trellis.par.get() # save them into a list called mysettings
names(mysettings)
# [1] "grid.pars"         "fontsize"          "background"        "panel.background" 
# [5] "clip"              "add.line"          "add.text"          "plot.polygon"     
# [9] "box.dot"           "box.rectangle"     "box.umbrella"      "dot.line"         
# [13] "dot.symbol"        "plot.line"         "plot.symbol"       "reference.line"   
# [17] "strip.background"  "strip.shingle"     "strip.border"      "superpose.line"   
# [21] "superpose.symbol"  "superpose.polygon" "regions"           "shade.colors"     
# [25] "axis.line"         "axis.text"         "axis.components"   "layout.heights"   
# [29] "layout.widths"     "box.3d"            "par.xlab.text"     "par.ylab.text"    
# [33] "par.zlab.text"     "par.main.text"     "par.sub.text"

mysettings$superpose.symbol
# $alpha
# [1] 1 1 1 1 1 1 1
# 
# $cex
# [1] 0.8 0.8 0.8 0.8 0.8 0.8 0.8
# 
# $col
# [1] "#0080ff"   "#ff00ff"   "darkgreen" "#ff0000"   "orange"    "#00ff00"   "brown"    
# 
# $fill
# [1] "#CCFFFF" "#FFCCFF" "#CCFFCC" "#FFE5CC" "#CCE6FF" "#FFFFCC" "#FFCCCC"
# 
# $font
# [1] 1 1 1 1 1 1 1
# 
# $pch
# [1] 1 1 1 1 1 1 1

# To change the default:
mysettings$superpose.symbol$pch <- c(1:10)
trellis.par.set(mysettings)

show.settings()


# 23.6. Customizing plot strips
library(lattice)
histogram(~height | voice.part, data = singer,
          strip = strip.custom(bg="lightgrey",
                               par.strip.text=list(col="black", cex=.8, font=3)),
          main = "Distribution of Heights by Voice Pitch",
          xlab = "Height (inches)")

# 23.7. Page arrangement
# figure 23.9
library(lattice)
graph1 <- histogram(~height | voice.part, data = singer,
                    main = "Heights of Choral Singers by Voice Part")
graph2 <- bwplot(height ~ voice.part, data = singer)
plot(graph1, split = c(1, 1, 1, 2))
plot(graph2, split = c(1, 2, 1, 2), newpage = FALSE)


# figure 23.10
library(lattice)
graph1 <- histogram(~height | voice.part, data = singer,
                    main = "Heights of Choral Singers by Voice Part")
graph2 <- bwplot(height ~ voice.part, data = singer)
plot(graph1, position = c(0, .33, 1, 1))
plot(graph2, position = c(0, 0, 1, .4), newpage = FALSE)


levels(singer$voice.part)
# [1] "Bass 2"    "Bass 1"    "Tenor 2"   "Tenor 1"   "Alto 2"    "Alto 1"    "Soprano 2"
# [8] "Soprano 1"

histogram(~height | voice.part, data = singer,
          index.cond=list(c(2,4,6,8,1,3,5,7)))
```
