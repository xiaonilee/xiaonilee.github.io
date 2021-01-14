rm(list = ls())
options(stringsAsFactors = F)

# load data
# setwd("/Users/xiaonili/Workspace/xiaonilee.github.io/content/post/frontoncol/fig8")
kegg <- read.csv("kegg2.csv", sep = ",", header = TRUE)
head(kegg)

mf <- read.csv("mf2.csv", sep = ",", header = TRUE)
head(mf)

bp <- read.csv("bp2.csv", sep = ",", header = TRUE)
head(bp)

cc <- read.csv("cc2.csv", sep = ",", header = TRUE)
head(cc)

library(ggplot2)

# Dot plot
p01 <- ggplot(kegg,aes(Count,reorder(Term,Count)))+
  geom_point(aes(size=Count,color=-1*log10(PValue))) +
  scale_color_gradient(low="red",high = "blue") +
  labs(color=expression(-log[10](PValue)), size="Count",
       x="Gene count",
       y=NULL,
       title = "KEGG")
ggsave("keggdotplot.png")

p02 <- ggplot(mf,aes(Count,reorder(Term,Count)))+
  geom_point(aes(size=Count,color=-1*log10(PValue))) +
  scale_color_gradient(low="red",high = "blue") +
  labs(color=expression(-log[10](PValue)), size="Count",
       x="Gene count",
       y="MF",
       title = "MF")
ggsave("mfdotplot.png")

p03 <- ggplot(cc,aes(Count,reorder(Term,Count)))+
  geom_point(aes(size=Count,color=-1*log10(PValue))) +
  scale_color_gradient(low="red",high = "blue") +
  labs(color=expression(-log[10](PValue)), size="Count",
       x="Gene count",
       y="CC",
       title = "CC")
ggsave("ccdotplot.png")

p04 <- ggplot(bp,aes(Count,reorder(Term,Count)))+
  geom_point(aes(size=Count,color=-1*log10(PValue))) +
  scale_color_gradient(low="red",high = "blue") +
  labs(color=expression(-log[10](PValue)), size="Count",
       x="Gene count",
       y="BP",
       title = "BP")
ggsave("bpdotplot.png")

require(gridExtra)
grid.arrange(p04, p03, p02, p01, ncol=1)


# Bar plot
p1 <- ggplot(data=kegg, aes(x=Count, y=reorder(Term,Count), fill=-log10(PValue))) +
  geom_bar(stat="identity",binwidth = 1) +
  scale_fill_gradient(low="red",high = "blue") +
  labs(color=expression(-1*log[10](PValue)),
       x="Gene count",
       y=NULL,
       title = "KEGG")
ggsave("keggbarplot.png")

p2 <- ggplot(data=mf, aes(x=Count, y=reorder(Term,Count), fill=-log10(PValue))) +
  geom_bar(stat="identity",binwidth = 1) +
  scale_fill_gradient(low="red",high = "blue") +
  labs(color=expression(-1*log[10](PValue)),
       x="Gene count",
       y="MF",
       title = "MF")
ggsave("mfbarplot.png")

p3 <- ggplot(data=cc, aes(x=Count, y=reorder(Term,Count), fill=-log10(PValue))) +
  geom_bar(stat="identity",binwidth = 1) +
  scale_fill_gradient(low="red",high = "blue") +
  labs(color=expression(-1*log[10](PValue)),
       x="Gene count",
       y="CC",
       title = "CC")
ggsave("ccbarplot.png")

p4 <- ggplot(data=bp, aes(x=Count, y=reorder(Term,Count), fill=-log10(PValue))) +
  geom_bar(stat="identity",binwidth = 1) +
  scale_fill_gradient(low="red",high = "blue") +
  labs(color=expression(-1*log[10](PValue)),
       x="Gene count",
       y="BP",
       title = "BP")
ggsave("bpbarplot.png")

require(gridExtra)
grid.arrange(p4, p3, p2, p1, ncol=1)
