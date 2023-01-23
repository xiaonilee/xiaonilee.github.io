
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

# Figure 8.1. 
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

# Figure 8.2
plot(women$height, women$weight,
     xlab = "Height (in inches)",
     ylab = "Weight (in pounds)")
lines(women$height, women$weight)


# with library(car) in the "car" package
# Figure 8.3.
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
# Figure 8.4.
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
# Figure 8.5.
plot(effect("hp:wt", fit,
            vcov. = vcov,
            list(wt=c(2.2,3.2,4.2))),
     multiline=TRUE)


# 8.3. Regression diagnostics
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data=states)
confint(fit)
# 2.5 %       97.5 %
#   (Intercept) -6.552191e+00 9.0213182149
# Population   4.136397e-05 0.0004059867
# Illiteracy   2.381799e+00 5.9038743192
# Income      -1.312611e-03 0.0014414600
# Frost       -1.966781e-02 0.0208304170


# 8.3.1. A typical approach
fit <- lm(weight ~ height, data=women)
par(mfrow=c(2,2))
# Figure 8.6.
plot(fit)

fit2 <- lm(weight ~ height + I(height^2), data=women)
par(mfrow=c(2,2))
# Figure 8.7.
plot(fit2)

fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data=states)
par(mfrow=c(2,2))
# Figure 8.8.
plot(fit)

# 8.3.2. An enhanced approach

# Normality
library(car)
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data=states)
# Figure 8.9.
qqPlot(fit, labels=row.names(states), id.method="identify",
       simulate=TRUE, main="Q-Q Plot")

states["Nevada",]
#         Murder Population Illiteracy Income Frost
# Nevada   11.5        590        0.5   5149   188

fitted(fit)["Nevada"]
#   Nevada 
# 3.878958 

residuals(fit)["Nevada"]
#   Nevada 
# 7.621042 

rstudent(fit)["Nevada"]
#   Nevada 
# 3.542929 

# Normality
# code listing 8.6. Function for plotting studentized residuals
residplot <- function(fit, nbreaks=10) {
  z <- rstudent(fit)
  hist(z, breaks=nbreaks, freq=FALSE,
       xlab="Studentized Residual",
       main="Distribution of Errors")
  rug(jitter(z), col="brown")
  curve(dnorm(x, mean=mean(z), sd=sd(z)),
        add=TRUE, col="blue", lwd=2)
  lines(density(z)$x, density(z)$y,
        col="red", lwd=2, lty=2)
  legend("topright",
         legend = c( "Normal Curve", "Kernel Density Curve"),
         lty=1:2, col=c("blue","red"), cex=.7)
}
# Figure 8.10.
residplot(fit)

# Independence of Errors
durbinWatsonTest(fit)
# lag Autocorrelation D-W Statistic p-value
# 1      -0.2006929      2.317691   0.306
# Alternative hypothesis: rho != 0


# Linearity
library(car)
# Figure 8.11.
crPlots(fit)

# Homoscedasticity
# code listing  8.7. Assessing homoscedasticity
library(car)
ncvTest(fit)
# Non-constant Variance Score Test 
# Variance formula: ~ fitted.values 
# Chisquare = 1.746514, Df = 1, p = 0.18632

# Figure 8.12.
spreadLevelPlot(fit)

# Suggested power transformation:  1.209626 


# 8.3.3. Global validation of linear model assumption
# code listing 8.8. Global test of linear model assumptions
# install.packages(gvlma)
library(gvlma)
gvmodel <- gvlma(fit)
summary(gvmodel)
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
# 
# 
# ASSESSMENT OF THE LINEAR MODEL ASSUMPTIONS
# USING THE GLOBAL TEST ON 4 DEGREES-OF-FREEDOM:
#   Level of Significance =  0.05 
# 
# Call:
#   gvlma(x = fit) 
# 
# Value p-value                Decision
# Global Stat        2.7728  0.5965 Assumptions acceptable.
# Skewness           1.5374  0.2150 Assumptions acceptable.
# Kurtosis           0.6376  0.4246 Assumptions acceptable.
# Link Function      0.1154  0.7341 Assumptions acceptable.
# Heteroscedasticity 0.4824  0.4873 Assumptions acceptable.


#============================================================

# 8.3.4. Multicollinearity
# code listing 8.9. Evaluating multicollinearity
library(car)
vif(fit)
# Population Illiteracy     Income      Frost 
#   1.245282   2.165848   1.345822   2.082547 


sqrt(vif(fit)) > 2
# Population Illiteracy     Income      Frost 
#      FALSE      FALSE      FALSE      FALSE 

#============================================================

# 8.4. Unusual observations

# 8.4.1. Outliers
library(car)
outlierTest(fit)

# 8.4.2. High leverage points

# Figure 8.13.
hat.plot <- function(fit) {
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  plot(hatvalues(fit), main="Index Plot of Hat Values")
  abline(h=c(2,3)*p/n, col="red", lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}
hat.plot(fit)

# 8.4.3. Influential observations

# Figure 8.14.
cutoff <- 4/(nrow(states)-length(fit$coefficients)-2)
plot(fit, which=4, cook.levels=cutoff)
abline(h=cutoff, lty=2, col="red")

# Figure 8.15. 
library(car)
avPlots(fit, ask=FALSE, onepage=TRUE, id.method="identify")

# Figure 8.16. Influence plot.
library(car)
influencePlot(fit, id.method="identify", main="Influence Plot",
              sub="Circle size is proportional to Cook's distance")


# 8.5.2. Transforming variables
# code listing 8.10. Box–Cox transformation to normality
library(car)
summary(powerTransform(states$Murder))
# bcPower Transformation to Normality 
# Est Power Rounded Pwr Wald Lwr Bnd Wald Upr Bnd
# states$Murder    0.6055           1       0.0884       1.1227
# 
# Likelihood ratio test that transformation parameter is equal to 0
# (log transformation)
# LRT df     pval
# LR test, lambda = (0) 5.665991  1 0.017297
# 
# Likelihood ratio test that no transformation is needed
# LRT df    pval
# LR test, lambda = (1) 2.122763  1 0.14512


boxTidwell(Murder~Population+Illiteracy,data=states)
#             MLE of lambda Score Statistic (z) Pr(>|z|)
# Population       0.86939             -0.3228   0.7468
# Illiteracy       1.35812              0.6194   0.5357
# 
# iterations =  19 

# 8.6.1. Comparing models
# code listing 8.11. Comparing nested models using the anova() function
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost,
           data=states)
fit2 <- lm(Murder ~ Population + Illiteracy, data=states)

anova(fit1,fit2)
# Analysis of Variance Table
# 
# Model 1: Murder ~ Population + Illiteracy + Income + Frost
# Model 2: Murder ~ Population + Illiteracy
# Res.Df    RSS Df Sum of Sq      F Pr(>F)
# 1     45 289.17                           
# 2     47 289.25 -2 -0.078505 0.0061 0.9939

# code listing 8.12. Comparing models with the AIC
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost,
           data=states)
fit2 <- lm(Murder ~ Population + Illiteracy, data=states)

AIC(fit1, fit2)
#       df      AIC
# fit1  6 241.6429
# fit2  4 237.6565


# 8.6.2. Variable selection
# code listing 8.13. Backward stepwise selection
library(MASS)
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost,
           data=states)
stepAIC(fit, direction = "backward")
# Start:  AIC=97.75
# Murder ~ Population + Illiteracy + Income + Frost
# 
# Df Sum of Sq    RSS     AIC
# - Frost       1     0.021 289.19  95.753
# - Income      1     0.057 289.22  95.759
# <none>                    289.17  97.749
# - Population  1    39.238 328.41 102.111
# - Illiteracy  1   144.264 433.43 115.986
# 
# Step:  AIC=95.75
# Murder ~ Population + Illiteracy + Income
# 
# Df Sum of Sq    RSS     AIC
# - Income      1     0.057 289.25  93.763
# <none>                    289.19  95.753
# - Population  1    43.658 332.85 100.783
# - Illiteracy  1   236.196 525.38 123.605
# 
# Step:  AIC=93.76
# Murder ~ Population + Illiteracy
# 
# Df Sum of Sq    RSS     AIC
# <none>                    289.25  93.763
# - Population  1    48.517 337.76  99.516
# - Illiteracy  1   299.646 588.89 127.311
# 
# Call:
#   lm(formula = Murder ~ Population + Illiteracy, data = states)
# 
# Coefficients:
#   (Intercept)   Population   Illiteracy  
# 1.6515497    0.0002242    4.0807366  


# code listing 8.14. All subsets regression
# install.packages("leaps")
library(leaps)
leaps <-regsubsets(Murder ~ Population + Illiteracy + Income +
                     Frost, data=states, nbest=4)
plot(leaps, scale="adjr2") # Figure 8.17

library(car)
subsets(leaps, statistic="cp",
        main="Cp Plot for All Subsets Regression")
abline(1,1,lty=2,col="red")


# 8.7.1. Cross-validation
# code listing 8.15. Function for k-fold cross-validated R-square
shrinkage <- function(fit, k=10){
  require(bootstrap)
  
  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
  
  x <- fit$model[,2:ncol(fit$model)]
  y <- fit$model[,1]
  
  results <- crossval(x, y, theta.fit, theta.predict, ngroup=k)
  r2 <- cor(y, fit$fitted.values)^2
  r2cv <- cor(y, results$cv.fit)^2
  cat("Original R-square =", r2, "\n")
  cat(k, "Fold Cross-Validated R-square =", r2cv, "\n")
  cat("Change =", r2-r2cv, "\n")
}

fit <- lm(Murder ~ Population + Income + Illiteracy + Frost, data=states)
shrinkage(fit)
# Original R-square = 0.5669502 
# 10 Fold Cross-Validated R-square = 0.4240783 
# Change = 0.1428719 


fit2 <- lm(Murder ~ Population + Illiteracy, data=states)
shrinkage(fit2)
# Original R-square = 0.5668327 
# 10 Fold Cross-Validated R-square = 0.5332059 
# Change = 0.03362679 


# 8.7.2. Relative importance
zstates <- as.data.frame(scale(states))
zfit <- lm(Murder~Population + Income + Illiteracy + Frost, data=zstates)
coef(zfit)
#   (Intercept)    Population        Income    Illiteracy         Frost 
# -2.054026e-16  2.705095e-01  1.072372e-02  6.840496e-01  8.185407e-03 


# code listing 8.16. relweights() function 
# for calculating relative importance of predictors
relweights <- function(fit,...){
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  barplot(t(import),names.arg=lbls,
          ylab="% of R-Square",
          xlab="Predictor Variables",
          main="Relative Importance of Predictor Variables",
          sub=paste("R-Square=", round(rsquare, digits=3)),
          ...)
  return(import)
}


relweights2 <- function(fit,...){
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  dotplot(t(import),names.arg=lbls,
          ylab="% of R-Square",
          xlab="Predictor Variables",
          main="Relative Importance of Predictor Variables",
          sub=paste("R-Square=", round(rsquare, digits=3)),
          ...)
  return(import)
}
# code listing 8.17. Applying the relweights() function
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data=states)
relweights(fit, col="lightgrey")
#             Weights
# Population 14.723401
# Illiteracy 59.000195
# Income      5.488962
# Frost      20.787442











































