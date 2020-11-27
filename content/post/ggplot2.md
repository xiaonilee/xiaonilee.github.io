---
title: "R Notebook : ggplot2"
date: 2020-09-19
lastmod: 2020-09-24
draft: false
tags: ["R", "Bioinformatics", "Plots", "R Package"]
categories: ["R", "Bioinformatics", "R Package"]
author: "Xiaoni"

weight: 1

mathjax: true

menu:
  main:
    parent: "docs"
    weight: 1
---

R package: ggplot2.

Data from Hadle Wickham 'ggplot2 Elegant Graphics for Data Analysis'.

<!--more-->

### Related Concept

```r
Data
Aesthetics
Geometries #几何对象
Facets
Statistics
Coordinates #坐标系
Themes
```

### Prerequisites

```r
install.packages("ggplot2")

library(ggplot2)
```

### Dataset

```r
data(mpg)
# if data is from outside: data <- read.table(file, header=T, row.names=NULL, sep=",)
data(economics)
```
  
  ![dataset](dataset.png)

### Basic Plot

#### 1. scatter

```r
ggplot (mpg, aes (x = displ, y = hwy)) +
geom_point()
```

![scatter](scatter.png)

#### 2. scatter add line

```r
ggplot(mpg, aes(displ, hwy)) +
geom_point() +
geom_smooth(method= "lm")
```

![scatter2](scatter2.png)

#### 3. scatter add color

```r
ggplot(mpg, aes (displ, hwy, colour=class)) +
geom_point()
```

![scatter3](scatter3.png)

#### 4. scatter and wrap
  
```r
ggplot(mpg, aes(displ, hwy)) +
geom_point() +
facet_wrap(~class)
```

![scatter4](scatter4.png)

#### 5. boxplot

```r
ggplot(mpg, aes(drv, hwy)) +
geom_boxplot()
```

![boxplot](boxplot.png)

#### 6. violin plot
  
```r
ggplot(mpg, aes(drv, hwy)) +
geom_violin()
```

![violin](violin.png)

#### 7. overlap
  
```r
ggplot(mpg, aes(drv, hwy)) +
geom_violin(trim=FALSE, fill='orange', color='red') +
geom_boxplot(width=0.1) + theme_minimal()
```

![overlap](overlap.png)

#### 8. histogram
  
```r
ggplot(mpg, aes(hwy)) +
geom_histogram()
```

![histogram](histogram.png)

#### 9. freqpoly
  
```r
ggplot(mpg, aes(hwy)) +
geom_freqpoly()
```

![freqpoly](freqpoly.png)

#### 10. bar
  
```r
ggplot(mpg, aes(manufacturer)) +
geom_bar()
```

![barplot](barplot.png)

#### 11. timeplot
  
```r
ggplot(economics, aes(date, unemploy / pop)) +
geom_line()
```

![timeplot](timeplot.png)

```r
ggplot(economics, aes(date, uempmed)) +
geom_line()
```

![timeplot2](timeplot2.png)

#### 12. path
  
```r
ggplot(economics, aes(unemploy / pop, uempmed)) +
geom_path() +
geom_point()
```

![path](path.png)

#### 13. contour
  
```r
ggplot(faithfuld, aes(eruptions, waiting)) +
geom_contour(aes(z = density, colour = ..level..))
```

![contour](contour.png)

#### 14. raster

```r
ggplot(faithfuld, aes(eruptions, waiting)) +
geom_raster(aes(fill = density))
```

![raster](raster.png)

### Fix with alpha

```r
norm <- ggplot(df, aes(x, y)) + xlab(NULL) + ylab(NULL)
norm + geom_point (alpha = 1 / 3)
norm + geom_point (alpha = 1 / 5)
norm + geom_point (alpha = 1 / 10)
```

### Setting theme

```r
different theme
theme_classic()
theme_minimal()
theme_economist()+ scale_colour_economist()
```

### ggplot2 extensions

- [ggthemes](https://www.rdocumentation.org/packages/ggthemes) Extra Themes, Scales and Geoms for 'ggplot2'.
- [ggpubr](https://www.rdocumentation.org/packages/ggpubr) Based Publication Ready Plots.
- [patchwork](https://www.rdocumentation.org/packages/patchwork) Expands the API, provide mathematical operators for combining multiple plots.
- [ggridges](https://www.rdocumentation.org/packages/ggridges) Provide a convenient way of visualizing changes in distributions over time or space.
- [ggdendro](https://www.rdocumentation.org/packages/ggdendro) Create Dendrograms and Tree Diagrams Using 'ggplot2'.
- [ggmap](https://www.rdocumentation.org/packages/ggmap) Spatial Visualization with ggplot2.
- [ggcorrplot](https://www.rdocumentation.org/packages/ggcorrplot) Visualize easily a correlation matrix using 'ggplot2'.
- [ggiraph](https://www.rdocumentation.org/packages/ggiraph) Create interactive 'ggplot2' graphics using 'htmlwidgets'.
- [GGally](https://www.rdocumentation.org/packages/GGaly) Extension to 'ggplot2',by adding several functions to reduce the complexity of combining geometric objects with transformed data.
- [ggalt](https://www.rdocumentation.org/packages/ggalt) Extra Coordinate Systems, Geoms and Statistical Transformations for 'ggplot2'.
- [ggforce](https://www.rdocumentation.org/packages/ggforce) Accelerating 'ggplot2'.
- [ggrepel](https://www.rdocumentation.org/packages/ggrepel) Automatically Position Non-Overlapping Text Labels with 'ggplot2'.
- [gganimate](https://www.rdocumentation.org/packages/gganimate) A Grammar of Animated Graphics.
- [ggradar](https://www.rdocumentation.org/packages/ggradar) Create radar charts.
- ggraph 绘制网络状、树状等特定形状的图形
- ggstance 实现常见图形的横向版本
- ggpmisc 光生物学相关扩展
- ggbio 提供基因组学数据图形
- geomnet 绘制网络状图形
- ggExtra 绘制图形的边界直方图
- plotROC 绘制交互式ROC曲线图
- ggspectra 绘制光谱图
- ggnetwork 网络状图形的geoms
- ggTimeSeries 时间序列数据可视化
- ggtree 树图可视化
- ggtern 绘制三元图
- ggseas 季节调整工具
- ggenealogy 浏览和展示系谱学数据

**In summary**, more knowledge by systematically learning [R related course](https://github.com/xiaonilee/Data_Analysis_with_R_byFacebook_ud651).
