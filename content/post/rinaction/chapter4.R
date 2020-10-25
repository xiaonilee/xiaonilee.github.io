
# Chapter4. Basic Data Management


# This chapter covers
#   Manipulating dates and missing values
#   Understanding data type conversions
#   Creating and recoding variables
#   Sorting, merging, and subsetting datasets
#   Selecting and dropping variables


# Remove most objects from the working environment
rm(list = ls()) 
options(stringsAsFactors = F)

# 4.1. A working example
# code listing 4.1. Creating the leadership data frame
manager <- c(1, 2, 3, 4, 5)
date <- c("10/24/08", "10/28/08", "10/1/08", "10/12/08", "5/1/09")
country <- c("US", "US", "UK", "UK", "UK")
gender <- c("M", "F", "F", "M", "F")
age <- c(32, 45, 25, 39, 99)
q1 <- c(5, 3, 3, 3, 2)
q2 <- c(4, 5, 5, 3, 2)
q3 <- c(5, 2, 5, 4, 1)
q4 <- c(5, 5, 5, NA, 2)
q5 <- c(5, 5, 2, NA, 1)
leadership <- data.frame(manager, date, country, gender, age,
                         q1, q2, q3, q4, q5, stringsAsFactors=FALSE)


# 4.2. Creating new variables
# code listing
mydata <- data.frame(x1 = c(2, 2, 6, 4),
                   x2 = c(3, 4, 2, 8))

mydata$sumx  <-  mydata$x1 + mydata$x2
mydata$meanx <- (mydata$x1 + mydata$x2)/2

attach(mydata)
mydata$sumx  <-  x1 + x2
mydata$meanx <- (x1 + x2)/2
detach(mydata)

mydata <- transform(mydata,
                    sumx  =  x1 + x2,
                    meanx = (x1 + x2)/2)
head(mydata)

# 4.3. Recoding variables
within()

# code
leadership <- within(leadership,{
  agecat <- NA
  agecat[age > 75]              <- "Elder"
  agecat[age >= 55 & age <= 75] <- "Middle Aged"
  agecat[age < 55]              <- "Young" })
# 4.4. Renaming variables
fix()

library(plyr)
rename()

# 4.5. Missing values
# NA, not available

is.na(leadership[,6:10])

# 4.5.2. Excluding missing values from analyses
x <- c(1, 2, NA, 3)
y <- sum(x, na.rm=TRUE)
y

# Deleting all observations with missing data (called listwise deletion) 
na.omit()

newdata <- na.omit(leadership)
head(newdata) # row 4th containing missing data is deleted.

# 4.6. Date values
as.Date()

mydates <- as.Date(c("2007-06-22", "2004-02-13"))
mydates

str(leadership)
myformat <- "%m/%d/%y"
leadership$date <- as.Date(leadership$date, myformat)
leadership$date
head(leadership)

Sys.Date()
# [1] "2020-10-22"

date()
# [1] "Thu Oct 22 09:46:54 2020"

today <- Sys.Date()
format(today, format="%B %d %Y")
# [1] "October 22 2020"
format(today, format="%A")
# [1] "Thursday"

today <- Sys.Date()
dob <- as.Date("2015-07-09")
difftime(today, dob, units = "weeks")

# Test the day in the week, same with today
difftime(today, dob, units = "days" )/7
# Time difference of 276 days

format(today, format="%A")
# [1] "Thursday"

# 4.7. Type conversions
# Test
is.numeric()
is.character()	
is.vector()	
is.matrix()	
is.data.frame()	
is.factor()	
is.logical()	

# Convert
as.numeric()
as.character()
as.vector()
as.matrix()
as.data.frame()
as.factor()
as.logical()

# 4.8. Sorting data
order()

# 4.9. Merging datasets

# 4.9.1. Adding columns horizontally
merge()
# by one or more common key variables
total <- merge(dataframeA, dataframeB, by="ID")

# without common key variables
total <- cbind(A, B)

# 4.9.2. Adding rows vertically
total <- rbind(dataframeA, dataframeB)

# 4.10. Subsetting datasets
# 4.10.1. Selecting (keeping) variables
myvars <- paste("q", 1:5, sep = "")
myvars
head(myvars)

# 4.10.2. Excluding (dropping) variables
# method1
myvars <- names(leadership) %in% c("q3", "q4")
newdata <- leadership[!myvars]

# equal method2
newdata <- leadership[c(-8,-9)]

# alternative method
leadership$q3 <- leadership$q4 <- NULL

# 4.10.3. Selecting observations

# code listing 4.6
newdata <- leadership[1:3,]

newdata <- leadership[which(leadership$gender=="M" &
                              leadership$age > 30),]

attach(leadership)
newdata <- leadership[which(gender=='M' & age > 30),]
detach(leadership)

# 4.10.4. The subset() function

# 4.10.5. Random samples
mysample <- leadership[sample(1:nrow(leadership), 3, replace=FALSE),]

# 4.11. Using SQL statements to manipulate data frames
# code listing 4.7 
library(sqldf)
newdf <- sqldf("select * from mtcars where carb=1 order by mpg",
               row.names=TRUE)
head(newdf)


sqldf("select avg(mpg) as avg_mpg, avg(disp) as avg_disp, gear 
      from mtcars where cyl in (4, 6) group by gear")
