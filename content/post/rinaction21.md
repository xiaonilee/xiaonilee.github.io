---
title: "Chapter 21. Creating a package"
date: 2020-11-20
lastmod: 2020-11-20
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

Notebook of Reading Books: R in Action_Chapter 21.

<!--more-->

## This chapter covers

- Creating the functions for a package

- Adding package documentation

- Building the package and sharing it with others

### 21.1. Nonparametric analysis and the npar package

- Figure 21.1. Distribution of healthy life expectancies at age 65 for women in the United States (2007–2009). 
  - The scores are **negatively skewed** (fewer scores at the low end).

  ![fig211](fig211.png)

- Figure 21.2. Dot chart of healthy life expectancies by region. 
  - The variability of HLE estimates differs across the four regions (compare the Northeast with the South).

  ![fig212](fig212.png)

#### 21.1.1. Comparing groups with the npar package

- **oneway {npar}**

- Figure 21.3. Annotated box plots displaying group differences. 
  - The plot is annotated with the medians and sample sizes for each group. 
  - The dotted vertical line represents the overall median.

  ![fig213](fig213.png)

### 21.2. Developing the package

#### 21.2.1. Computing the statistics

#### 21.2.2. Printing the results

#### 21.2.3. Summarizing the results

#### 21.2.4. Plotting the results

#### 21.2.5. Adding sample data to the package

#### 21.3. Creating the package documentation

- Table 21.2. Tags for use with Roxygen2

  ![tab1](tab1.png)

### 21.4. Building the package

- 1. Install the necessary tools.

- 2.Set up the directories.
  - [oneway.R](oneway.R), [print.R](print.R), [summary.R](summary.R), [plot.R](plot.R), [npar.R](npar.R), and [life.R](life.R) are in the R file.
  ![rfile](rfil.png)
  
  - [life.rda](life.rda) is in the data file.
  ![dat](dat.png)

- 3. Generate the documentation.
  
  ![genrt](genrt.png)

  ![tree1](tree1.png)

  ![tree2](tree2.png)

- 4. Bulid the package.
  
  ![bui1](build1.png)

- 5. Check the package(optional).

  ![che](check1.png)
  ![che](check2.png)

- 6. Create a [PDF](npar.pdf) manual(optional).

- 7. Install the package locally (optional).

  ![inst](insta.png)

- 8. Upload the package to CRAN(optional).

Attach is the [Script](chapter21.R) of chapter21.

Show me the code <i class="far fa-hand-pointer"></i>

```r
# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)


pkg <- "npar_1.0.tar.gz"
loc <- "http://www.statmethods.net/RiA"
url <- paste(loc, pkg, sep="/")
download.file(url, pkg)
install.packages(pkg, repos=NULL, type="source")


# 21.1. Nonparametric analysis and the npar package
library(npar)
hist(life$hlef, xlab = "Healthy Life Expectancy (years) at Age 65",
     main = "Distribution of Healthy Life Expectancy for Women",
     col = "grey", breaks = 10) # figure 21.1

# figure 21.2
library(ggplot2)
ggplot(data = life, aes(x=region, y=hlef)) +
  geom_point(size=3, color="darkgrey") +
  labs(title = "Distribution of HLE Estimates by Region",
       x="US Region", y="Healthy Life Expectancy at Age 65") +
  theme_bw()


# 21.1.1. Comparing groups with the npar package
# code listing 21.1. Comparison of HLE estimates with the npar package
library(npar)
results <- oneway(hlef ~ region, life)
summary(results)
# data: hlef on region 
# 
# Omnibus Test
# Kruskal-Wallis chi-squared = 17.8749, df = 3, p-value = 0.0004668
# 
# Descriptive Statistics
#          South North Central    West Northeast
# n      16.0000      12.00000 13.0000   9.00000
# median 13.0000      15.40000 15.6000  15.70000
# mad     1.4826       1.26021  0.7413   0.59304
# 
# Multiple Comparisons (Wilcoxon Rank Sum Tests)
# Probability Adjustment = holm
#         Group.1       Group.2    W           p   
# 1         South North Central 28.0 0.008583179 **
# 2         South          West 27.0 0.004737844 **
# 3         South     Northeast 17.0 0.008583179 **
# 4 North Central          West 63.5 1.000000000   
# 5 North Central     Northeast 42.0 1.000000000   
# 6          West     Northeast 54.5 1.000000000   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


plot(results, col="lightblue", main="Multiple Comparisons",
     xlab="US Region",
     ylab="Healthy Life Expectancy (years) at Age 65")


#==========================================================
# 21.2.1. Computing the statistics
# code listing 21.2. Contents of the oneway.R file

#' @title Nonparametric group comparisons
#'
#' @description
#' \code{oneway} computes nonparametric group comparisons, including an
#' omnibus test and post-hoc pairwise group comparisons.
#'
#' @details
#' This function computes an omnibus Kruskal-Wallis test that the
#' groups are equal, followed by all pairwise comparisons using
#' Wilcoxon Rank Sum tests. Exact Wilcoxon tests can be requested if
#' there are no ties on the dependent variable. The p-values are
#' adjusted for multiple comparisons using the \code{\link{p.adjust}}
#' function.
#'
#' @param formula an object of class formula, relating the dependent
#' variable to the grouping variable.
#' @param data a data frame containing the variables in the model.
#' @param exact logical. If \code{TRUE}, calculate exact Wilcoxon tests.
#' @param sort logical. If \code{TRUE}, sort groups by median dependent
#' variable values.
#' @param method method for correcting p-values for multiple comparisons.
#' @export
#' @return a list with 7 elements:
#' \item{CALL}{function call}
#' \item{data}{data frame containing the depending and grouping variable}
#' \item{sumstats}{data frame with descriptive statistics by group}
#' \item{kw}{results of the Kruskal-Wallis test}
#' \item{method}{method used to adjust p-values}
#' \item{wmc}{data frame containing the multiple comparisons}
#' \item{vnames}{variable names}
#' @author Rob Kabacoff <rkabacoff@@statmethods.net>
#' @examples
#' results <- oneway(hlef ~ region, life)
#' summary(results)
#' plot(results, col="lightblue", main="Multiple Comparisons",
#'      xlab="US Region", ylab="Healthy Life Expectancy at Age 65")

# 1. Function call
oneway <- function(formula, data, exact=FALSE, sort=TRUE,	
                   method=c("holm", "hochberg", "hommel", "bonferroni",	
                            "BH", "BY", "fdr", "none")){	
  
  # 2. Checks arguments
  if (missing(formula) || class(formula) != "formula" ||	
      length(all.vars(formula)) != 2)	
    stop("'formula' is missing or incorrect")	
  
  method <- match.arg(method)	
  
  # 3. Sets up data
  df <- model.frame(formula, data)	
  y <- df[[1]]	
  g <- as.factor(df[[2]])	
  vnames <- names(df)	
  
  # 4. Reorders factor levels
  if(sort) g <- reorder(g, y, FUN=median)	
  groups <- levels(g)
  k <- nlevels(g)
  
  # 5. Summary statistics
  getstats <- function(x)(c(N = length(x), Median = median(x),	
                            MAD = mad(x)))	
  sumstats <- t(aggregate(y, by=list(g), FUN=getstats)[2])	
  rownames(sumstats) <- c("n", "median", "mad")	
  colnames(sumstats) <- groups	
  
  # 6. Statistical tests
  kw <- kruskal.test(formula, data)	
  wmc <- NULL	
  for (i in 1:(k-1)){	
    for (j in (i+1):k){	
      y1 <- y[g==groups[i]]	
      y2 <- y[g==groups[j]]	
      test <- wilcox.test(y1, y2, exact=exact)	
      r <- data.frame(Group.1=groups[i], Group.2=groups[j],	
                      W=test$statistic[[1]], p=test$p.value)	
      # note the [[]] to return a single number	
      wmc <- rbind(wmc, r)	
    }	
  }	
  wmc$p <- p.adjust(wmc$p, method=method)	
  
  # 7.Return results
  data <- data.frame(y, g)	
  names(data) <- vnames	
  results <- list(CALL = match.call(),	
                  data=data,	
                  sumstats=sumstats, kw=kw,	
                  method=method, wmc=wmc, vnames=vnames)	
  class(results) <- c("oneway", "list")	
  return(results)	
}


# 21.2.2. Printing the results
print(results)
# data: hlef by region 
# 
# Multiple Comparisons (Wilcoxon Rank Sum Tests)
# Probability Adjustment = holm
#         Group.1       Group.2    W           p
# 1         South North Central 28.0 0.008583179
# 2         South          West 27.0 0.004737844
# 3         South     Northeast 17.0 0.008583179
# 4 North Central          West 63.5 1.000000000
# 5 North Central     Northeast 42.0 1.000000000
# 6          West     Northeast 54.5 1.000000000


# code listing 21.3. Contents of the print.R file
#' @title Print multiple comparisons
#'
#' @description
#' \code{print.oneway} prints pairwise group comparisons.
#'
#' @details
#' This function prints Wilcoxon pairwise multiple comparisons created
#' by the \code{\link{oneway}} function.
#'
#' @param x an object of class \code{oneway}.
#' @param ... additional arguments passed to the function.
#' @method print oneway
#' @export
#' @return the input object is returned silently.
#' @author Rob Kabacoff <rkabacoff@@statmethods.net>
#' @examples
#' results <- oneway(hlef ~ region, life)
#' print(results)
print.oneway <- function(x, ...){
  # 1. Checks input
  if (!inherits(x, "oneway"))	
    stop("Object must be of class 'oneway'")	
  
  # 2.Print the header
  cat("data:", x$vnames[1], "by", x$vnames[2], "\n\n")	
  cat("Multiple Comparisons (Wilcoxon Rank Sum Tests)\n")	
  cat(paste("Probability Adjustment = ", x$method, "\n", sep=""))	
  
  # 3. print the table
  print(x$wmc,  ...)	
}



# 21.2.3. Summarizing the results
summary(results)

# code lisitng 21.4. Contents of the summary.R file
#' @title Summarize oneway nonparametric analyses
#'
#' @description
#' \code{summary.oneway} summarizes the results of a oneway
#' nonparametric analysis.
#'
#' @details
#' This function prints a summary of analyses produced by
#' the \code{\link{oneway}} function. This includes descriptive
#' statistics by group, an omnibus Kruskal-Wallis test, and
#' Wilcoxon pairwise multiple comparisons.
#'
#' @param object an object of class \code{oneway}.
#' @param ... additional parameters.
#' @method summary oneway
#' @export
#' @return the input object is returned silently.
#' @author Rob Kabacoff <rkabacoff@@statmethods.net>
#' @examples
#' results <- oneway(hlef ~ region, life)
#' summary(results)
summary.oneway <- function(object, ...){
  if (!inherits(object, "oneway"))
    stop("Object must be of class 'oneway'")
  
  if(!exists("digits")) digits <- 4L
  
  kw <- object$kw
  wmc <- object$wmc
  cat("data:", object$vnames[1], "on", object$vnames[2], "\n\n")
  
  
  # Kruskal-Wallis test
  cat("Omnibus Test\n")	
  cat(paste("Kruskal-Wallis chi-squared = ",	
            round(kw$statistic,4),	
            ", df = ", round(kw$parameter, 3),	
            ", p-value = ",	
            format.pval(kw$p.value, digits = digits),	
            "\n\n", sep=""))	
  
  # Descriptive Statistics
  cat("Descriptive Statistics\n")	
  print(object$sumstats, ...)	
  
  # Table annotation
  wmc$stars <- " "	
  wmc$stars[wmc$p <   .1] <- "."	
  wmc$stars[wmc$p <  .05] <- "*"	
  wmc$stars[wmc$p <  .01] <- "**"	
  wmc$stars[wmc$p < .001] <- "***"	
  names(wmc)[which(names(wmc)=="stars")] <- " "	
  
  # Pairwise multiple comparisons
  cat("\nMultiple Comparisons (Wilcoxon Rank Sum Tests)\n")	
  cat(paste("Probability Adjustment = ", object$method, "\n", sep=""))	
  print(wmc, ...)	
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' '	
     1\n")	
}


# 21.2.4. Plotting the results
# code listing 21.5. Contents of the plot.R file
#' @title Plot nonparametric group comparisons
#'
#' @description
#' \code{plot.oneway} plots nonparametric group comparisons.
#'
#' @details
#' This function plots nonparametric group comparisons
#' created by the \code{\link{oneway}} function using
#' annotated side by side boxplots. Medians and
#' sample sizes are placed at the top of the chart.
#' The overall median is represented by a horizontal
#' dashed line.
#'
#' @param x an object of class \code{oneway}.
#' @param ... additional arguments passed to the
#' \code{\link{boxplot}} function.
#' @method plot oneway
#' @export
#' @return NULL
#' @author Rob Kabacoff <rkabacoff@@statmethods.net>
#' @examples
#' results <- oneway(hlef ~ region, life)
#' plot(results, col="lightblue", main="Multiple Comparisons",
#'      xlab="US Region", ylab="Healthy Life Expectancy at Age 65")
plot.oneway <- function(x, ...){
  
  # 1. Checks input
  if (!inherits(x, "oneway"))	
    stop("Object must be of class 'oneway'")	
  
  # 2. Generates the box plots
  data <- x$data	
  y <- data[,1]	
  g <- data[,2]	
  stats <- x$sumstats	
  lbl <- paste("md=", stats[2,], "\nn=", stats[1,], sep="")	
  opar <- par(no.readonly=TRUE)	
  par(mar=c(5,4,8,2))
  boxplot(y~g,  ...)
  
  # 3. Annotates the plot
  abline(h=median(y), lty=2, col="darkgrey")	
  axis(3, at=1:length(lbl), labels=lbl, cex.axis=.9)	
  par(opar)
}


# 21.2.5. Adding sample data to the package
# code listing 21.6. Creating the life data frame
region <- c(rep("North Central", 12), rep("Northeast", 9),
            rep("South", 16), rep("West", 13))

state <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI",
           "CT","ME","MA","NH","NJ","NY","PA","RI","VT","AL","AR","DE",
           "FL","GA","KY","LA","MD","MS","NC","OK","SC","TN","TX","VA",
           "WV","AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT",
           "WA","WY")
hlem <- c(12.6,12.2,13.4,13.1,12.8,14.3,11.7,13.1,12.9,12.2,13.3,13.4,
          14.3,13.5,13.8,14,12.9,13.6,12.8,13.1,13.9,10.3,11.6,13.5,
          14.3,11.6,10.2,11.6,13.3,10.1,11.7,10.8,12,11.2,12.2,13.3,
          10.3,13.3,13.7,13.8,14.3,15,13.1,13.4,12.8,13.1,13.9,14.3,14,
          13.7)
hlef <- c(14.3,14.1,15.9,15.1,14.7,16.7,14,15.7,16,14,16.4,16.1,16.7,
          15.7,15.9,16,14.8,15.3,14.8,15.6,16.2,11.7,12.7,15.7,16.4,
          13.1,11.6,12.3,15.3,11.4,13.5,12.9,13.6,12.5,13.4,14.9,11.6,
          14.9,16.3,15.5,16.2,17.3,15.1,15.6,14.5,14.7,16,15.7,16,15.2)

life <- data.frame(region=factor(region), state=factor(state), hlem, hlef)

save(life, file='life.rda')


# code listing 21.7. Contents of the life.R file
#' @title Healthy Life Expectancy at Age 65
#'
#' @description A dataset containing the healthy life expectancy (expected
#' years of life in good health) at age 65, by US state in 2007-2009.
#' Estimates are reported separately for men and women.
#'
#' @docType data
#' @keywords datasets
#' @name life
#' @usage life
#' @format A data frame with 50 rows and 4 variables. The variables
#' are as follows:
#' \describe{
#'   \item{region}{A factor with 4 levels (North Central, Northeast,
#'                 South, West)}
#'   \item{state}{A factor with the 2-letter ISO codes for the 50 US
#'                states}
#'   \item{hlem}{Healthy life expectancy for men in years}
#'   \item{hlef}{Healthy life expectancy for women in years}
#' }
#' @source The \code{hlem} and \code{hlef} data were obtained from
#' the Center for Disease Control and Prevention
#' \emph{Morbidity and Mortality Weekly Report} at \url{
#' http://www.cdc.gov/mmwr/preview/mmwrhtml/mm6228a1.htm?s_cid=mm6228a1_w}.
#' The \code{region} variable was added from the
#' \code{\link[datasets]{state.region}} dataset.
NULL


# 21.3. Creating the package documentation
# code listing 21.8. Contents of the npar.R file
#' Functions for nonparametric group comparisons.
#'
#' npar provides tools for calculating and visualizing
#' nonparametric differences among groups.
#'
#' @docType package
#' @name npar-package
#' @aliases npar
NULL

#... this file must end with a blank line after the NULL...
```
