
# Chapter 22. Creating dynamic reports


# This chapter covers
#   Publishing results to the web
#   Incorporating R results into Microsoft Word or Open Document reports
#   Creating dynamic reports, where changing the data changes the report
#   Creating publication quality documents with R, Markdown, and LaTeX


# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)


# 22.2. Creating dynamic reports with R and Markdown
# code listing 22.1. The document of women.Rmd: a Markdown template with embedded R code

# render the file
# setwd() to current path
library(rmarkdown)
render("women.Rmd", "html_document")

# render the file
library(rmarkdown)
# A new version of TeX Live has been released, so
# tinytex::reinstall_tinytex()
render("womenPDF.Rmd", "pdf_document")

# render the file
library(rmarkdown)
render("womenWord.Rmd", "word_document")


# 22.3. Creating dynamic reports with R and LaTeX
# code listing 22.2. drugs.Rnw: a sample LaTeX template with embedded R code

library(knitr)
knit("drugs.Rnw")

library(knitr)
knit2pdf("drugs.Rnw")