library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)


myRDS <- readRDS("15days_preprocessed.rds") 

groups <- sample(c("15dayS1", "15dayS2"), size =14364 ,replace = TRUE)
names(groups) <- colnames(myObject)
myObject  <-  AddMetaData(object = myObject, metadata = groups, col.name = "sample")
saveRDS(myObject, file = "15daySplit_preprocessed.rds")

#After analysing and merging 
myRDS <- readRDS("photoreceptor_analysed.rds")
groups <- sample(c("15dayS1", "15dayS2"), size =27095 ,replace = TRUE)
groups <- sample(c("15dayS1", "15dayS2", "30dayS1", "30dayS2"), size =27095 ,replace = TRUE)
names(groups) <- colnames(myObject)
myObject  <-  AddMetaData(object = myObject, metadata = groups, col.name = "sample")
saveRDS(myObject, file = "photoreceptor_split.rds")





