---
title: "Chapter 19. Advanced graphics with ggplot2"
date: 2020-11-17
lastmod: 2020-11-17
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

Notebook of Reading Books: R in Action_Chapter 19.

<!--more-->

## This chapter covers

- An introduction to the ggplot2 package

- Using shape, color, and size to visualize multivariate data

- Comparing groups with faceted graphs

- Customizing ggplot2 plots

### 19.1. The four graphics system in R

  ![tab1](tab1.png)

### 19.2. An introduction to ggplot2 package

- Figure 19.1. Scatterplot of automobile weight by mileage.

  ![fig191](fig191.png)

- Figure 19.2. Scatterplot of automobile weight by gas mileage, with a superimposed line of best fit and 95% confidence region.
  
  ![fig192](fig192.png)

- Figure 19.3. A scatterplot showing the relationship between horsepower and gas mileage separately for transmission and engine type. 
  - The number of cylinders in each automobile engine is represented by both shape and color.

  ![fig193](fig193.png)

### 19.3. Specifying the plot type with geoms

- Table 19.2. Geom functions
  
  ![tab2](tab2.png)

- Table 19.3. Common options for geom functions

  ![tab3](tab3.png)

- Figure 19.4. Histogram of singer heights.

  ![fig194](fig194.png)

- Figure 19.5. Box plot of singer heights by voice part.
  
  ![fig195](fig195.png)

- Figure 19.6. Notched box plots with superimposed points describing the salaries of college professors by rank. 
  - A rug plot is provided on the vertical axis.

  ![fig196](fig196.png)

- Figure 19.7. A combined violin and box plot graph of singer heights by voice part.

  ![fig197](fig197.png)

### 19.4. Grouping

- Figure 19.8. Density plots of university salaries, grouped by academic rank.
  
  ![fig198](fig198.png)

- Figure 19.9. Scatterplot of years since graduation and salary. 
  - Academic rank is represented by color, 
  - and sex is represented by shape.

  ![fig199](fig199.png)

- Figure 19.10. Three versions of a grouped bar chart. 
  - Each displays the number of professors by academic rank and sex.

![fig1910](fig1910.png)

- Figure 19.11. Faceted graph showing the distribution (histogram) of singer heights by voice part.

  ![fig1911](fig1911.png)

- Figure 19.12. Scatterplot of years since graduation and salary. 
  - Academic rank is represented by color and shape, and sex is faceted.

  ![fig1912](fig1912.png)

- Figure 19.13. Faceted density plots for singer heights by voice part. 

  ![fig1913](fig1913.png)

### 19.6. Adding smoothed lines

- Figure 19.14. Scatterplot of years since doctorate and current faculty salary. 
  - A fitted loess smoothed line with 95% confidence limits has been added.

  ![fig1914](fig1914.png)

- Figure 19.15. Scatterplot of years since graduation vs. salary with separate fitted quadratic regression lines for men and women.

  ![fig1915](fig1915.png)

### 19.7. Modifying the appearance of ggplot2 graphs

#### 19.7.1. Axes

  ![tab4](tab4.png)

- Figure 19.16. Box plots of faculty salaries grouped by academic rank and sex. 
  - The axis text has been customized.

  ![fig1916](fig1916.png)

#### 19.7.2. Legends

- Figure 19.17. Box plots of faculty salaries grouped by academic rank. 
  - The axis text has been customized, along with the ***legend title and position***.

  ![fig1917](fig1917.png)

#### 19.7.3 Scales

- Figure 19.18. Bubble chart of auto weight by mileage, with point size representing engine displacement.

  ![fig1918](fig1918.png)

- Figure 19.19. Scatterplot of salary vs. experience for assistant, associate, and full professors.
  - Point colors have been specified ***manually***.

  ![fig1919](fig1919.png)

- Figure 19.19-2. Scatterplot of salary vs. experience for assistant, associate, and full professors.
  - Point colors have been specified ***manually*** with `scale_color_brewer()` function.

- Figure 19.19-3. colors sets with `display.brewer.all()`

  ![fig1919-3](fig1919-3.png)

#### 19.7.4. Themes

- Figure 19.20. Box plots with **a customized theme**. 

  ![fig1920](fig1920.png)

#### 19.7.5. Multiple graphs per page

- Figure 19.21. Placing three ggplot2 plots in a single graph with `grid.arrange()`.

  ![fig1921](fig1921.png)

#### Saving graphs

- `ggsave()`

Attach is the [Script](chapter19.R) of chapter19.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)


# 19.2. An introduction to the ggplot2 package
# figure 19.1
library(ggplot2)
ggplot(data = mtcars, aes(x=wt, y=mpg)) +
  geom_point() +
  labs(title = "Automobile Data", x="Weight", y="Miles Per Gallon")


# figure 19.2
library(ggplot2)
ggplot(data=mtcars, aes(x=wt, y=mpg)) +
  geom_point(pch=17, color="blue", size=2) +
  geom_smooth(method="lm", color="red", linetype=2) +
  labs(title="Automobile Data", x="Weight", y="Miles Per Gallon")

mtcars$am <- factor(mtcars$am, levels=c(0,1),
                    labels=c("Automatic", "Manual"))
mtcars$vs <- factor(mtcars$vs, levels=c(0,1),
                    labels=c("V-Engine", "Straight Engine"))
mtcars$cyl <- factor(mtcars$cyl)


# figure 19.3
library(ggplot2)
ggplot(data=mtcars, aes(x=hp, y=mpg,
                        shape=cyl, color=cyl)) +
  geom_point(size=3) +
  facet_grid(am~vs) +
  labs(title="Automobile Data by Engine Type",
       x="Horsepower", y="Miles Per Gallon")


# 19.3. Specifying the plot type with geoms
# figure 19.4
data(singer, package = "lattice")
ggplot(singer, aes(x=height)) + geom_histogram()


# figure 19.5
ggplot(singer, aes(x=voice.part, y=height)) + geom_boxplot()


# figure 19.6
data(Salaries, package="carData")

dim(Salaries)
# [1] 397   6

head(Salaries)

#        rank discipline yrs.since.phd yrs.service  sex salary
# 1      Prof          B            19          18 Male 139750
# 2      Prof          B            20          16 Male 173200
# 3  AsstProf          B             4           3 Male  79750
# 4      Prof          B            45          39 Male 115000
# 5      Prof          B            40          41 Male 141500
# 6 AssocProf          B             6           6 Male  97000

colnames(Salaries)
# [1] "rank"          "discipline"    "yrs.since.phd" "yrs.service"   "sex"          
# [6] "salary" 
library(ggplot2)
ggplot(Salaries, aes(x=rank, y=salary)) +
  geom_boxplot(fill="cornflowerblue",
               color="black", notch=TRUE) +
  geom_point(position="jitter", color="red", alpha=.5) +
  geom_rug(side="l", color="black")


# figure 19.7
library(ggplot2)
data(singer, package="lattice")
ggplot(singer, aes(x=voice.part, y=height)) +
  geom_violin(fill="lightpink") +
  geom_boxplot(fill="lightgray", width=.2)


# 19.4. Grouping
# figure 19.8
data(Salaries, package="carData")
library(ggplot2)
ggplot(data=Salaries, aes(x=salary, fill=rank)) +
  geom_density(alpha=.3)

# figure 19.9
ggplot(Salaries, aes(x=yrs.since.phd, y=salary, color=rank,
                     shape=sex)) + geom_point()


# figure 19.10
require(gridExtra)
plot1 <- ggplot(Salaries, aes(x=rank, fill=sex)) +
  geom_bar(position="stack") + labs(title='position="stack"')

plot2 <- ggplot(Salaries, aes(x=rank, fill=sex)) +
  geom_bar(position="dodge") + labs(title='position="dodge"')

plot3 <- ggplot(Salaries, aes(x=rank, fill=sex)) +
  geom_bar(position="fill") + labs(title='position="fill"')

grid.arrange(plot1, plot2, plot3, ncol=3)


# 19.5. Faceting
# figure 19.11
data(singer, package="lattice")
library(ggplot2)
ggplot(data=singer, aes(x=height)) +
  geom_histogram() +
  facet_wrap(~voice.part, nrow=4)


# figure 19.12
library(ggplot2)
ggplot(Salaries, aes(x=yrs.since.phd, y=salary, color=rank,
                     shape=rank)) + geom_point() + facet_grid(.~sex)


# figure 19.13
data(singer, package="lattice")
library(ggplot2)
ggplot(data=singer, aes(x=height, fill=voice.part)) +
  geom_density() +
  facet_grid(voice.part~.)


# 19.6. Adding smoothed lines
# figure 19.14
data(Salaries, package="carData")
library(ggplot2)
ggplot(data=Salaries, aes(x=yrs.since.phd, y=salary)) +
  geom_smooth() +
  geom_point()

# figure 19.15
ggplot(data=Salaries, aes(x=yrs.since.phd, y=salary,
                          linetype=sex, shape=sex, color=sex)) +
  geom_smooth(method=lm, formula=y~poly(x,2),
              se=FALSE, size=1) +
  geom_point(size=2)

# 19.7.1. Axes
# figure 19.16
data(Salaries,package="carData")
library(ggplot2)
ggplot(data=Salaries, aes(x=rank, y=salary, fill=sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("AsstProf", "AssocProf", "Prof"),
                   labels=c("Assistant\nProfessor",
                            "Associate\nProfessor",
                            "Full\nProfessor")) +
  scale_y_continuous(breaks=c(50000, 100000, 150000, 200000),
                     labels=c("$50K", "$100K", "$150K", "$200K")) +
  labs(title="Faculty Salary by Rank and Sex", x="", y="")

# 19.7.2. Legends
# figure 19.17
data(Salaries,package="carData")
library(ggplot2)
ggplot(data=Salaries, aes(x=rank, y=salary, fill=sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks=c("AsstProf", "AssocProf", "Prof"),
                   labels=c("Assistant\nProfessor",
                            "Associate\nProfessor",
                            "Full\nProfessor")) +
  scale_y_continuous(breaks=c(50000, 100000, 150000, 200000),
                     labels=c("$50K", "$100K", "$150K", "$200K")) +
  labs(title="Faculty Salary by Rank and Gender",
       x="", y="", fill="Gender") +
  theme(legend.position=c(.1,.8))

# 19.7.3 Scales
# figure 19.18
ggplot(mtcars, aes(x=wt, y=mpg, size=disp)) +
  geom_point(shape=21, color="black", fill="cornsilk") +
  labs(x="Weight", y="Miles Per Gallon",
       title="Bubble Chart", size="Engine\nDisplacement")


# figure 19.19
data(Salaries, package="carData")
ggplot(data=Salaries, aes(x=yrs.since.phd, y=salary, color=rank)) +
  scale_color_manual(values=c("orange", "olivedrab", "navy")) +
  geom_point(size=2)

# figure 19.19-2
ggplot(data=Salaries, aes(x=yrs.since.phd, y=salary, color=rank)) +
  scale_color_brewer(palette="Set1") + 
  geom_point(size=2)

library(RColorBrewer)
display.brewer.all()


# 19.7.4. Themes
data(Salaries, package="carData")
library(ggplot2)
mytheme <- theme(plot.title=element_text(face="bold.italic",
                                         size="14", color="brown"),
                 axis.title=element_text(face="bold.italic",
                                         size=10, color="brown"),
                 axis.text=element_text(face="bold", size=9,
                                        color="darkblue"),
                 panel.background=element_rect(fill="white",
                                               color="darkblue"),
                 panel.grid.major.y=element_line(color="grey",
                                                 linetype=1),
                 panel.grid.minor.y=element_line(color="grey",
                                                 linetype=2),
                 panel.grid.minor.x=element_blank(),
                 legend.position="top")

ggplot(Salaries, aes(x=rank, y=salary, fill=sex)) +
  geom_boxplot() +
  labs(title="Salary by Rank and Sex", x="Rank", y="Salary") +
  mytheme

# alternative no define of mytheme
ggplot(Salaries, aes(x=rank, y=salary, fill=sex)) +
  geom_boxplot() +
  labs(title="Salary by Rank and Sex", x="Rank", y="Salary") +
  theme(plot.title=element_text(face="bold.italic",
                                size="14", color="brown"),
        axis.title=element_text(face="bold.italic",
                                size=10, color="brown"),
        axis.text=element_text(face="bold", size=9,
                               color="darkblue"),
        panel.background=element_rect(fill="white",
                                      color="darkblue"),
        panel.grid.major.y=element_line(color="grey",
                                        linetype=1),
        panel.grid.minor.y=element_line(color="grey",
                                        linetype=2),
        panel.grid.minor.x=element_blank(),
        legend.position="top")


# 19.7.5. Multiple graphs per page
data(Salaries, package="carData")
library(ggplot2)
p1 <- ggplot(data=Salaries, aes(x=rank)) + geom_bar()
p2 <- ggplot(data=Salaries, aes(x=sex)) + geom_bar()
p3 <- ggplot(data=Salaries, aes(x=yrs.since.phd, y=salary)) + geom_point()

library(gridExtra)
grid.arrange(p1, p2, p3, ncol=3)

# Saving graphs
myplot <- ggplot(data=mtcars, aes(x=mpg)) + geom_histogram()
ggsave(file="mygraph.png", plot=myplot, width=5, height=4)
```
