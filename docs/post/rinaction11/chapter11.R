
# Chapter 11. Intermediate graphs


# This chapter covers
#   Visualizing bivariate and multivariate relationships
#   Working with scatter and line plots
#   Understanding correlograms
#   Using mosaic and association plots


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













