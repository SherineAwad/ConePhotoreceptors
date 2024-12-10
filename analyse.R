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
myRDS <- args[1]

split_string <- strsplit(myRDS, ".rds")[[1]]
mysample <- split_string[1]
print(mysample)

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"


all.genes <- rownames(myObject)
myObject <- NormalizeData(myObject)
myObject <- FindVariableFeatures(myObject)
myObject<- ScaleData(myObject, features = all.genes)
myObject <- RunPCA(myObject,reduction.key = "PC_")
myObject <- FindNeighbors(myObject, dims = 1:8)
myObject <- FindClusters(myObject, resolution = 2.0)
myObject <- RunUMAP( myObject, dims = 1:8, reduction.name = "umap") 
myObject[["RNA"]] <- JoinLayers(myObject[["RNA"]])

myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)

