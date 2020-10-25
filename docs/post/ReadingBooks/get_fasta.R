rm(list = ls())
options(stringsAsFactors = F)

library(GEOquery)
gpl21827 <- getGEO("GPL21827", destdir = '.')
save(gpl21827, file = 'GPL21827.Rdata')

load("GPL21827.Rdata")

class(gpl21827)
length(gpl21827)

gpl21827
colnames(Table(gpl21827))
probe2seq <- Table(gpl21827)[,c(1,4)]
head(probe2seq)

# x1 = probe2seq[1,]
# paste0('>', x1[1], '\n')

res <- paste(apply(probe2seq, 1, function(x) 
  paste0('>', x[1],'\n', x[2])), collapse = '\n')
head(res)

temp <- tempfile()
tempfile()

write(res, temp)



