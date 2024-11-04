library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myRDS <- paste(mysample, "_analysed.rds", sep="")
myRDS

myObject <- readRDS(myRDS)


figure_name <- ""
figure_name <- paste(mysample, "_UMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "sample",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

gene_list <- rownames(myObject)



