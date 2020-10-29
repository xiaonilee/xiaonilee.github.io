
# Chapter 8. Regression


# This chapter covers
#   Fitting and interpreting linear models
#   Evaluating model assumptions
#   Selecting among competing models


# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)


# 8.2.1. Fitting regression models with lm()
myfit <- lm(formula, data)

# 8.2.2. Simple linear regression
# code listing Simple linear regression
fit <- lm(weight ~ height, data = women)
colnames(women)
summary(fit)
# Call:
#   lm(formula = weight ~ height, data = women)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.7333 -1.1333 -0.3833  0.7417  3.1167 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -87.51667    5.93694  -14.74 1.71e-09 ***
#   height        3.45000    0.09114   37.85 1.09e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.525 on 13 degrees of freedom
# Multiple R-squared:  0.991,	Adjusted R-squared:  0.9903 
# F-statistic:  1433 on 1 and 13 DF,  p-value: 1.091e-14


fitted(fit)
women$weight
residuals(fit)

plot(women$height, women$weight,
     xlab = "Height (in inches)",
     ylab = "Weight (in pounds)")
abline(fit)

# 8.2.3. Polynomial regression
fit2 <- lm(weight ~ height + I(height^2), data=women)

# code listing 8.2. Polynomial regression
fit2 <- lm(weight ~ height + I(height^2), data = women)
summary(fit2)
# Call:
#   lm(formula = weight ~ height + I(height^2), data = women)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.50941 -0.29611 -0.00941  0.28615  0.59706 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 261.87818   25.19677  10.393 2.36e-07 ***
#   height       -7.34832    0.77769  -9.449 6.58e-07 ***
#   I(height^2)   0.08306    0.00598  13.891 9.32e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3841 on 12 degrees of freedom
# Multiple R-squared:  0.9995,	Adjusted R-squared:  0.9994 
# F-statistic: 1.139e+04 on 2 and 12 DF,  p-value: < 2.2e-16


plot(women$height, women$weight,
     xlab = "Height (in inches)",
     ylab = "Weight (in pounds)")
lines(women$height, women$weight)


# with library(car) in the "car" package
library(car)
scatterplot(weight ~ height,
            data=women,
            spread=FALSE, smoother=list(lty=2),
            pch=19,
            main="Women Age 30-39",
            xlab="Height (inches)",
            ylab="Weight (lbs.)")

# 8.2.4. Multiple linear regression

# code listing 8.3 Examining bivariate relationships
states <- as.data.frame(state.x77[,c("Murder", "Population",
                                     "Illiteracy", "Income", "Frost")])

cor(states)
library(car)
scatterplotMatrix(states, spread=F, lty.smooth=2,
                  main = "Scatter Plot Matrix")


# code listing 8.4 Multiple linear regression
states <- as.data.frame(state.x77[,c("Murder", "Population",
                                     "Illiteracy", "Income", "Frost")])

fit <- lm(Murder ~ Population + Illiteracy + Income + Frost,
          data = states)
summary(fit)
# Call:
#   lm(formula = Murder ~ Population + Illiteracy + Income + Frost, 
#      data = states)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.7960 -1.6495 -0.0811  1.4815  7.6210 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.235e+00  3.866e+00   0.319   0.7510    
# Population  2.237e-04  9.052e-05   2.471   0.0173 *  
#   Illiteracy  4.143e+00  8.744e-01   4.738 2.19e-05 ***
#   Income      6.442e-05  6.837e-04   0.094   0.9253    
# Frost       5.813e-04  1.005e-02   0.058   0.9541    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.535 on 45 degrees of freedom
# Multiple R-squared:  0.567,	Adjusted R-squared:  0.5285 
# F-statistic: 14.73 on 4 and 45 DF,  p-value: 9.133e-08

#===================================================================
# 8.2.5. Multiple linear regression with interactions
# code listing 8.5 Multiple linear regression with a significant interaction term
fit <- lm(mpg ~ hp + wt + hp:wt, data=mtcars)
summary(fit)
# Call:
#   lm(formula = mpg ~ hp + wt + hp:wt, data = mtcars)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.0632 -1.6491 -0.7362  1.4211  4.5513 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 49.80842    3.60516  13.816 5.01e-14 ***
#   hp          -0.12010    0.02470  -4.863 4.04e-05 ***
#   wt          -8.21662    1.26971  -6.471 5.20e-07 ***
#   hp:wt        0.02785    0.00742   3.753 0.000811 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.153 on 28 degrees of freedom
# Multiple R-squared:  0.8848,	Adjusted R-squared:  0.8724 
# F-statistic: 71.66 on 3 and 28 DF,  p-value: 2.981e-13

#====================================================================

# install.packages("effects")
library(effects)
plot(effect("hp:wt", fit,
            vcov. = vcov,
            list(wt=c(2.2,3.2,4.2))),
     multiline=TRUE)

























