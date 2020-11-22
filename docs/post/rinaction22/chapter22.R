
# Chapter 22. Creating dynamic reports


# This chapter covers
#   Publishing results to the web
#   Incorporating R results into Microsoft Word or Open Document reports
#   Creating dynamic reports, where changing the data changes the report
#   Creating publication quality documents with R, Markdown, and LaTeX


rm(list = ls())
options(stringsAsFactors = F)


# 22.2. Creating dynamic reports with R and Markdown
# code listing 22.1. women.Rmd: a Markdown template with embedded R code

library(rmarkdown)
render("/Users/xiaonili/Workspace/xiaonilee.github.io/content/post/rinaction22/women.Rmd", "html_document")


library(rmarkdown)
# A new version of TeX Live has been released, so
# tinytex::reinstall_tinytex()
render("/Users/xiaonili/Workspace/xiaonilee.github.io/content/post/rinaction22/womenPDF.Rmd", "pdf_document")











