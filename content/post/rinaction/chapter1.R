# 10/19/2020 Day1
# Chapter1. Introduction To R


# This chapter covers
#   Installing R
#   Understanding the R language
#   Running programs


# Remove most objects from the working environment
rm(list = ls()) 
options(stringsAsFactors = F)

# Section 1.3.2 Getting help
help.start()

## Section 1.3.3 Functions for managing the R workspace
getwd()
setwd("mydirectory")
ls()            # List the objects in the current workspace.
rm(objectlist)  # Remove (delete) one or more objects.
options()       # View or set current options.
history()       # Display your last # commands (default = 25).
save()
load()
q()

## Section 1.3.4 Input and Output

# Input
source("filename")

source("script1.R") # Working through an example.

# Output
sink("myoutput", append=TRUE, split=TRUE)
pdf("mygraphs.pdf")
source("script2.R")

sink()
dev.off()
source("script3.R")