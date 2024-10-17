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

mysample = args[1] 
nRNA1 = as.double(args [2]) 
nRNA2 = as.double(args[3])
features = as.double(args[4])
mt = as.double(args [5])

nRNA1
nRNA2
features 
mt

myRDS = paste(mysample, "_preprocessed.rds", sep="") 
myRDS
myObject <- readRDS(myRDS)


myObject <- subset(x = myObject,subset = nCount_RNA < nRNA1 & nCount_RNA > nRNA2 & nFeature_RNA > features & percent.mt < mt)

myRDS <- paste(mysample, "_filtered.rds", sep="") 
saveRDS(myObject, file = myRDS)



