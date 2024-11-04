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


myObject <- readRDS("15days_filtered.rds") 

myObject@meta.data[,'sample']<-apply(as.matrix(rownames(myObject@meta.data)), 1, function(X1){X2<-strsplit(X1,'-')[[1]][2];}) 
myObject@meta.data[,'sample'] <- recode(myObject@meta.data[,'sample'], "1" = "15dayS1", "2" = "15dayS2")
saveRDS(myObject, file = "15days_split.rds")


myObject1 <- readRDS("30dayS1_filtered.rds")
myObject2 <- readRDS("30dayS2_filtered.rds")

myObject <- merge(myObject1, y = c(myObject2) )

myObject@meta.data[,'sample']<-apply(as.matrix(rownames(myObject@meta.data)), 1, function(X1){X2<-strsplit(X1,'-')[[1]][2];})
myObject@meta.data[,'sample'] <- recode(myObject@meta.data[,'sample'], "1_1" = "30dayS1", "1_2" = "30dayS2")

saveRDS(myObject, file = "30days_split.rds")



