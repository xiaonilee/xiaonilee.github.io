
# Chapter 18. Advanced methods for missing data


# This chapter covers
#   Identifying missing data
#   Visualizing missing data patterns
#   Complete-case analysis
#   Multiple imputation of missing data

# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 18.2. Identifying missing values
# install.packages("VIM")
data(sleep, package="VIM")

# list the rows that do not have missing values
sleep[complete.cases(sleep),]

# list the rows that have one or more missing values
sleep[!complete.cases(sleep),]


sum(is.na(sleep$Dream))
# [1] 12

mean(is.na(sleep$Dream))
# [1] 0.1935484

mean(!complete.cases(sleep))
# [1] 0.3225806

# 18.3.1. Tabulating missing values
# install.packages("mice")
library(mice)
data(sleep, package = "VIM")
md.pattern(sleep) # figure 18.2

# 18.3.2. Exploring missing data visually
library("VIM")
aggr(sleep, prop=FALSE, numbers=TRUE) # figure 18.2-2

matrixplot(sleep, sortby = "BrainWgt") # figure 18.3

marginplot(sleep[c("Gest", "Dream")], pch = c(20),
           col = c("darkgray", "red", "blue")) # figure 18.4

# 18.3.3. Using correlations to explore missing values
x <- as.data.frame(abs(is.na(sleep)))

head(sleep, 5)
#    BodyWgt BrainWgt NonD Dream Sleep Span Gest Pred Exp Danger
# 1 6654.000   5712.0   NA    NA   3.3 38.6  645    3   5      3
# 2    1.000      6.6  6.3   2.0   8.3  4.5   42    3   1      3
# 3    3.385     44.5   NA    NA  12.5 14.0   60    1   1      1
# 4    0.920      5.7   NA    NA  16.5   NA   25    5   2      3
# 5 2547.000   4603.0  2.1   1.8   3.9 69.0  624    3   5      4

head(x, 5)
#   BodyWgt BrainWgt NonD Dream Sleep Span Gest Pred Exp Danger
# 1       0        0    1     1     0    0    0    0   0      0
# 2       0        0    0     0     0    0    0    0   0      0
# 3       0        0    1     1     0    0    0    0   0      0
# 4       0        0    1     1     0    1    0    0   0      0
# 5       0        0    0     0     0    0    0    0   0      0

y <- x[which(apply(x, 2, sum) > 0)]

cor(y)
#              NonD       Dream       Sleep        Span        Gest
# NonD   1.00000000  0.90711474  0.48626454  0.01519577 -0.14182716
# Dream  0.90711474  1.00000000  0.20370138  0.03752394 -0.12865350
# Sleep  0.48626454  0.20370138  1.00000000 -0.06896552 -0.06896552
# Span   0.01519577  0.03752394 -0.06896552  1.00000000  0.19827586
# Gest  -0.14182716 -0.12865350 -0.06896552  0.19827586  1.00000000

cor(sleep, y, use = "pairwise.complete.obs")
#                 NonD       Dream        Sleep        Span        Gest
# BodyWgt   0.22682614  0.22259108  0.001684992 -0.05831706 -0.05396818
# BrainWgt  0.17945923  0.16321105  0.007859438 -0.07921370 -0.07332961
# NonD              NA          NA           NA -0.04314514 -0.04553485
# Dream    -0.18895206          NA -0.188952059  0.11699247  0.22774685
# Sleep    -0.08023157 -0.08023157           NA  0.09638044  0.03976464
# Span      0.08336361  0.05981377  0.005238852          NA -0.06527277
# Gest      0.20239201  0.05140232  0.159701523 -0.17495305          NA
# Pred      0.04758438 -0.06834378  0.202462711  0.02313860 -0.20101655
# Exp       0.24546836  0.12740768  0.260772984 -0.19291879 -0.19291879
# Danger    0.06528387 -0.06724755  0.208883617 -0.06666498 -0.20443928


#=======================================================================
# 18.6. Complete-case analysis(listwise deletion)
options(digits = 1)
cor(na.omit(sleep))
#          BodyWgt BrainWgt NonD Dream Sleep  Span  Gest  Pred  Exp Danger
# BodyWgt     1.00     0.96 -0.4 -0.07  -0.3  0.47  0.71  0.10  0.4   0.26
# BrainWgt    0.96     1.00 -0.4 -0.07  -0.3  0.63  0.73 -0.02  0.3   0.15
# NonD       -0.39    -0.39  1.0  0.52   1.0 -0.37 -0.61 -0.35 -0.6  -0.53
# Dream      -0.07    -0.07  0.5  1.00   0.7 -0.27 -0.41 -0.40 -0.5  -0.57
# Sleep      -0.34    -0.34  1.0  0.72   1.0 -0.38 -0.61 -0.40 -0.6  -0.60
# Span        0.47     0.63 -0.4 -0.27  -0.4  1.00  0.65 -0.17  0.3   0.01
# Gest        0.71     0.73 -0.6 -0.41  -0.6  0.65  1.00  0.09  0.6   0.31
# Pred        0.10    -0.02 -0.4 -0.40  -0.4 -0.17  0.09  1.00  0.6   0.93
# Exp         0.41     0.32 -0.6 -0.50  -0.6  0.32  0.57  0.63  1.0   0.79
# Danger      0.26     0.15 -0.5 -0.57  -0.6  0.01  0.31  0.93  0.8   1.00

# alternative
# cor(sleep, use="complete.obs)


#========================================================================
fit <- lm(Dream ~ Span + Gest, data = na.omit(sleep))
summary(fit)
# Call:
# lm(formula = Dream ~ Span + Gest, data = na.omit(sleep))
# 
# Residuals:
#    Min     1Q Median     3Q    Max 
# -2.333 -0.915 -0.221  0.382  4.183 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.480122   0.298476    8.31  3.7e-10 ***
# Span        -0.000472   0.013130   -0.04    0.971    
# Gest        -0.004394   0.002081   -2.11    0.041 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1 on 39 degrees of freedom
# Multiple R-squared:  0.167,	Adjusted R-squared:  0.125 
# F-statistic: 3.92 on 2 and 39 DF,  p-value: 0.0282


# 18.7. Multiple imputation
library(mice)
data(sleep, package = "VIM")
imp <- mice(sleep, seed = 1234)

fit <- with(imp, lm(Dream ~ Span + Gest))
pooled <- pool(fit)
summary(pooled)
#          term estimate std.error statistic df p.value
# 1 (Intercept)    2.597     0.249      10.4 52   2e-14
# 2        Span   -0.004     0.012      -0.3 56   7e-01
# 3        Gest   -0.004     0.001      -3.0 55   5e-03



imp
# Class: mids
# Number of multiple imputations:  5 
# Imputation methods:
#   BodyWgt BrainWgt     NonD    Dream    Sleep     Span     Gest     Pred      Exp   Danger 
#        ""       ""    "pmm"    "pmm"    "pmm"    "pmm"    "pmm"       ""       ""       "" 
# PredictorMatrix:
#          BodyWgt BrainWgt NonD Dream Sleep Span Gest Pred Exp Danger
# BodyWgt        0        1    1     1     1    1    1    1   1      1
# BrainWgt       1        0    1     1     1    1    1    1   1      1
# NonD           1        1    0     1     1    1    1    1   1      1
# Dream          1        1    1     0     1    1    1    1   1      1
# Sleep          1        1    1     1     0    1    1    1   1      1
# Span           1        1    1     1     1    0    1    1   1      1
# Number of logged events:  5 
#   it im  dep meth   out
# 1  3  2 Span  pmm Sleep
# 2  3  2 Gest  pmm Sleep
# 3  4  2 Span  pmm Sleep
# 4  4  2 Gest  pmm Sleep
# 5  4  4 Span  pmm Sleep

imp$imp$Dream
#      1   2   3   4   5
# 1  0.0 0.5 0.5 0.5 0.3
# 3  0.5 1.4 1.5 1.5 1.3
# 4  3.6 4.1 3.1 4.1 2.7
# 14 0.3 1.0 0.5 0.0 0.0
# 24 3.6 0.8 1.4 1.4 0.9
# 26 2.4 0.5 3.9 3.4 1.2
# 30 2.6 0.8 2.4 2.2 3.1
# 31 0.6 1.3 1.2 1.8 2.1
# 47 1.3 1.8 1.8 1.8 3.9
# 53 0.5 0.5 0.6 0.5 0.3
# 55 2.6 3.6 2.4 1.8 0.5
# 62 1.5 3.4 3.9 3.4 2.2












