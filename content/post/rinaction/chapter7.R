
# Chapter7. Basic statistics


# This chapter covers
#   Descriptive statistics
#   Frequency and contingency tables
#   Correlations and covariances
#   t-tests
#   Nonparametric statistics


# Remove most objects from the working environment
rm(list = ls()) 
options(stringsAsFactors = F)


# 7.1. Descriptive statisticss

vars <- c("mpg", "hp", "wt")
head(mtcars[vars])

# 7.1.1. A menagerie of methods

# code listing 7.1  Descriptive statistics via summary()
summary(mtcars[vars])

# code listing 7.2 Descriptive statistics via sapply()
mystats <- function(x, na.omit=F) {
              if (na.omit)
                x <- x[!is.na(x)]
              m <- mean(x)
              n <- length(x)
              s <- sd(x)
              skew <- sum((x-m)^3/s^3)/n
              kurt <- sum((x-m)^4/s^4)/n -3
              return(c(n=n, mean=m, stdev=s, skew=skew, kurtosis=kurt))
}

sapply(mtcars[vars], mystats)
sapply(mtcars[vars], mystats, na.omit=T)

# 7.1.1.1 Extensions
# code listing 7.3 Descriptive statistics via describe() in the Hmisc package()
library(Hmisc)
describe(mtcars[vars])

# code listing 7.4 Descriptive statistics via stat.desc() in the pastecs package
# install.packages("pastecs")
library(pastecs)
stat.desc(mtcars[vars])

# code listing 7.5 Descriptive statistics via describe() in the psych package
# install.packages("psych")
library(psych)
describe(mtcars[vars])

# 7.1.2 Descriptive statistics by group
# code listing 7.6 Descriptive statistics by group using aggregate()
aggregate(mtcars[vars], by=list(am=mtcars$am), mean)
mtcars$am

aggregate(mtcars[vars], by=list(am=mtcars$am), sd)

# code listing 7.7 Descriptive statistics by group using by()
dstats <- function(x) sapply(x, mystats)
by(mtcars[vars], mtcars$am, dstats)

# 7.1.2.1 Extensions
# code listing 7.8 Summary statistics by group using summaryBy() in the doBy package
# install.packages("doBy")
library(doBy)
summaryBy(mpg+hp+wt~am, data=mtcars, FUN=mystats)

# code listing 7.9 Summary statistics by group using describe.by() in the psych package
library(psych)
describeBy(mtcars[vars], mtcars$am)
describeBy(mtcars[vars], list(am=mtcars$am))

# 7.1.3 Visualizing results

# 7.2. Frequency and contingency tables
library(vcd)

# 7.2.1 Generating frequency tables

# One-Way Tables
mytable <- with(Arthritis,table(Improved))
mytable

prop.table(mytable)
prop.table(mytable)*100

# Two-Way Tables

mytable <- xtabs(~ Treatment+Improved, data = Arthritis)
mytable
table(Arthritis$Treatment)

margin.table(mytable,1)
margin.table(mytable,2)
prop.table(mytable,1)
prop.table(mytable,2)

addmargins(mytable)
addmargins(prop.table(mytable))
addmargins(prop.table(mytable, 1), 2)
addmargins(prop.table(mytable, 2), 1)

# code listing 7.10 Two-way table using CrossTable
# install.packages("gmodels")
library(gmodels)
CrossTable(Arthritis$Treatment, Arthritis$Improved)
?CrossTable

# Multidimensional Tables
# code listing Listing 7.11. Three-way contingency table
mytable <- xtabs(~ Treatment+Sex+Improved, data = Arthritis)
mytable
ftable(mytable)


margin.table(mytable, 1)
margin.table(mytable, 2)
margin.table(mytable, 3)

margin.table(mytable, c(1, 3))
margin.table(mytable, c(1, 2))
ftable(prop.table(mytable, c(1, 2)))
ftable(addmargins(prop.table(mytable, c(1, 2)), 3))


# 7.2.2. Tests of independence
# The significance tests in the section evaluated 
#       whether or not sufficient evidence existed 
#       to reject a null hypothesis of independence between variables

# Chi-Square Test of Independence
# code listing 7.13. Chi-square test of independence
library(vcd)
mytable <- xtabs(~Treatment+Improved, data = Arthritis)
chisq.test(mytable)

mytable <- xtabs(~Treatment+Sex, data = Arthritis)
chisq.test(mytable)

# Fisher’s Exact Test
mytable <- xtabs(~Treatment+Improved, data = Arthritis)
fisher.test(mytable)

# Cochran–Mantel–Haenszel Test
mytable <- xtabs(~Treatment+Sex+Improved, data = Arthritis)
mantelhaen.test(mytable)

# 7.2.3. Measures of association

# code Listing 7.14. Measures of association for a two-way table
library(vcd)
mytable <- xtabs(~Treatment+Improved, data = Arthritis)
assocstats(mytable)

# 7.2.4. Visualizing results

# 7.2.5. Converting tables to flat files
# code Listing 7.15. Converting a table into a flat file via table2flat
table2flat <- function(mytable) {
  df <- as.data.frame(mytable)
  rows  <- dim(df)[1]
  cols  <- dim(df)[2]
  x <- NULL
  for (i in 1:rows){
    for (j in 1:df$Freq[i]){
      row <- df[i,c(1:(cols-1))]
      x <- rbind(x,row)
    }
  }
  row.names(x)<-c(1:dim(x)[1])
  return(x)
  
}

# code Listing 7.16. Using the table2flat() function with published data
treatment <- rep(c("Placebo", "Treated"), times=3)
improved <- rep(c("None", "Some", "Marked"), each=2)
Freq <- c(29,13,7,17,7,21)
mytable <- as.data.frame(cbind(treatment, improved, Freq))
mydata <- table2flat(mytable)
head(mydata)

# 7.3. Correlations
# 7.3.1. Types of correlations
# Pearson, Spearman, and Kendall Correlations
# code listing 7.17. Covariances and correlations
states <- state.x77[, 1:6]
colnames(state.x77)
dim(state.x77)

cov(states)
cor(states)
cor(states, method = "spearman")
cor(states, method = "kendall")
cor(states, method = "pearson")

x <- states[,c("Population", "Income", "Illiteracy", "HS Grad")]
y <- states[,4:5]
cor(x,y)

# Partial Correlations
# install.packages("ggm")

library(ggm)
# partial correlation of population and murder rate, controlling
# for income, illiteracy rate, and HS graduation rate
pcor(c(1,5,2,3,6), cov(states))

# Other Types of Correlations


# 7.3.2. Testing correlations for significance

# to test an individual Pearson, Spearman, and Kendall correlation coefficient
cor.test()

# code listing Testing a correlation coefficient for significance
cor.test(states[,3], states[,5])

# code listing 7.19. Correlation matrix and tests of significance via corr.test
library(psych)
corr.test(states, use = "complete")

# 7.3.3. Visualizing correlations

# 7.4. t-tests

# 7.4.1. Independent t-test
library(MASS)
t.test(Prob ~ So, data = UScrime)


# 7.4.2. Dependent t-test
library(MASS)
colnames(UScrime)

sapply(UScrime[c("U1", "U2")], function(x) (c(mean=mean(x), sd=sd(x))))
with(UScrime, t.test(U1, U2, paired=T))

# 7.4.2. Dependent t-test

# 7.5. Nonparametric tests of group differences

# 7.5.1. Comparing two groups
with(UScrime, by(Prob, So, median))

wilcox.test(Prob ~ So, data=UScrime)

sapply(UScrime[c("U1", "U2")], median)
median(UScrime$U1)
median(UScrime$U2)

with(UScrime, wilcox.test(U1, U2, paired=TRUE))

# 7.5.2. Comparing more than two groups
# When there are more than two groups to be compared, 
# you must turn to other methods

#=======dataset,dataset=====useful method=======cbind()======#
states <- as.data.frame(cbind(state.region, state.x77))

kruskal.test(Illiteracy ~ state.region, data=states)

# 7.6. Visualizing group differences

# 7.7. Summary


















































