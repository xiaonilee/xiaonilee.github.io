
# Chapter2. Creating Dataset


# This chapter covers
#   Exploring R data structures
#   Using data entry
#   Importing data
#   Annotating datasets


# Remove most objects from the working environment
rm(list = ls()) 
options(stringsAsFactors = F)

help.start()

# Section 2.2 Data structures
attach(); detach(); with() 
# example
summary(mtcars$mpg)
plot(mtcars$mpg, mtcars$disp)
plot(mtcars$mpg, mtcars$wt)

# equal:
attach(mtcars)
summary(mpg)
plot(mpg, disp)
plot(mpg, wt)
detach(mtcars)

# equal:
with(mtcars, {
  summary(mpg, disp, wt)
  plot(mpg, disp)
  plot(mpg, wt)
})

## Section 2.2.6 Lists
# code listing 2.7 Creating a list
g <- "My First List"
h <- c(25, 26, 18, 39)
j <- matrix(1:10, nrow = 5)
k <- c("one", "two", "three")
mylist <- list(title=g, ages=h, j, k)
mylist

str(mylist)
mylist[[2]]
mylist[[1]]

# column
mylist[[3]][,2]

# row
mylist[[3]][2,]
mylist[[4]][1]

# A Note about 'block comments'
if(FALSE) {...}


## Section 2.3.1 Data Input

## Section 2.3.1 Entering Data from keyboard
mydatatxt <- "
age gender weight
25 m 166
30 f 115
18 f 120
"
mydata <- read.table(header = TRUE, text = mydatatxt)
head(mydata)
str(mydata)
class(mydata)
summary(mydata)

write.table(mydata,"mydata.csv", col.names = TRUE,
            row.names = TRUE, append = FALSE)

ls() # display current objects

# method 2 
mydata <- data.frame(age=numeric(0),
                     gender=character(0), weight=numeric(0))
mydata <- edit(mydata)


# 2.5. Useful functions for working with data objects

dim(object)	# Dimensions of an object.
str(object)	# Structure of an object.
class(object)	# Class or type of an object.
head(object)	# Lists the first part of the object.
tail(object)	# Lists the last part of the object.
c(object, object,...)	# Combines objects into a vector.
cbind(object, object, ...)	# Combines objects as columns.
rbind(object, object, ...)	# Combines objects as rows.

rm(object, object, ...)	# Deletes one or more objects. 
rm(list = Is()) # remove most objects from the working environment.