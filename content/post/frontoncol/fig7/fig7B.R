rm(list = ls())
options(stringsAsFactors = F)

# Reading the data
a <- read.table('/Users/xiaonili/Downloads/fig7222.txt',header = FALSE)
class(a)

# remove row1 and row2
mydata <- a[c(3:14),]

# add rowname
rownames(mydata) <- mydata$V1

# remove col1
mydata <- mydata[, c(2:539)]

# check colnames and rownames
colnames(mydata)
rownames(mydata)


# Transposing the dataframe
## method1. Converting it into matrix
mydata_transpose <- as.data.frame(t(as.matrix(mydata)))

## method2. Using the data.table library
library(data.table)
mydata_transpose2 <- transpose(mydata)
rownames(mydata_transpose2) <- colnames(mydata)
colnames(mydata_transpose2) <- rownames(mydata)

# data.frame 2 numeric
df <- as.data.frame(sapply(mydata_transpose, function(x) as.numeric(as.character(x))))
df2 <- as.matrix(df)

# # put together
# DF = as.matrix(as.data.frame(lapply(mydata_transpose, as.numeric)))

fig7B <- cor(df2, use = "complete.obs", method = "spearman")
write.csv()




