rm(list = ls())
options(stringsAsFactors = F)

## Load library
library(openxlsx)
library(tidyverse)

library(limma)
library(readr)

## load data
TCGA_rawdata <- read_tsv("/Users/xiaonili/Downloads/TCGA-HNSC.htseq_counts.tsv.gz")

dim(TCGA_rawdata)

probeMap <- read.table("/Users/xiaonili/Downloads/gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T)
probeMap[1:4,1:4]

## ID reverse
TCGA_gset <- TCGA_rawdata %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>%
  select(gene, starts_with("TCGA") )
TCGA_gset[1:4,1:4]


## average replicate genes
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene))

colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)
write.csv(TCGA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_Countdata_log2+1.csv")
TCGA_gset[1:4,1:4]

## group by patient.id
TCGA_group_list <- ifelse(as.numeric(substring(colnames(TCGA_gset),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)


# recognize mRNA lncRNA and miRNA

## load data
mRNA_info <- read.xlsx("/Users/xiaonili/Downloads/Gene_info.xlsx",sheet = "mRNA_info")
lncRNA_info <- read.xlsx("/Users/xiaonili/Downloads/Gene_info.xlsx",sheet = "lncRNA_info")
miRNA_info <- read.xlsx("/Users/xiaonili/Downloads/Gene_info.xlsx",sheet = "miRNA_info")

## Get data.matrix for mRNA
mRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% mRNA_info$gene_name,]
dim(mRNA_gset)


write.csv(mRNA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_mRNA.csv",quote = F,row.names = T)

## Get data.matrix for lncRNA
lncRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% lncRNA_info$gene_name,]
dim(lncRNA_gset)

write.csv(lncRNA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_lncRNA.csv",quote = F,row.names = T)

## Get data.matrix for miRNA
miRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% miRNA_info$gene_name,]
dim(miRNA_gset)


write.csv(miRNA_gset,"/Users/xiaonili/Downloads/TCGA_HNSC_miRNA.csv",quote = F,row.names = T)


# Match clinical and survival information with expression

## Load clinical data
Phenodata <- read_tsv("/Users/xiaonili/Downloads/TCGA-HNSC.GDC_phenotype.tsv.gz")

Phenodata[1:4,1:4]

Phenodata$submitter_id.samples <- substring(Phenodata$submitter_id.samples,1,15) %>% 
  gsub("-",".",.)
Phenodata[1:4,1:4]

##Load survival data
Sur_data <- read_tsv("/Users/xiaonili/Downloads/TCGA-HNSC.survival.tsv.gz")

Sur_data$sample <- substring(Sur_data$sample,1,15) %>% gsub("-",".",.)
Sur_data[1:4,1:4]

## merge data and choose interested col
Phen_surv <- Phenodata %>%
  inner_join(Sur_data,by = c("submitter_id.samples" = "sample")) %>%
  select(submitter_id.samples,age_at_index.demographic,gender.demographic,
         tumor_grade.diagnoses,neoplasm_histologic_grade,tumor_stage.diagnoses,OS,OS.time)
head(Phen_surv)

## match expression with phenodata and do order
Phen_surv = Phen_surv[match(colnames(TCGA_gset),Phen_surv$submitter_id.samples),]
identical(Phen_surv$submitter_id.samples,colnames(TCGA_gset))

## group
Phen_surv$group <- TCGA_group_list
Phen_surv = dplyr::select(Phen_surv,submitter_id.samples,group,everything())
write.csv(Phen_surv,"/Users/xiaonili/Downloads/TCGA_HNSC_phenotype.csv")
head(Phen_surv)




