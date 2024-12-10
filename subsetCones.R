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
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
myRDS <- args[1]

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]
print(mysample)

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"

cones <- subset(myObject, idents = "Cones", invert = FALSE)
table(cones@active.ident)
table(cones@meta.data[,'sample'])
cones@active.ident <- as.factor(cones@meta.data[,'sample'])
names(cones@active.ident) <- rownames(cones@meta.data)

cones <- FindNeighbors(cones, dims = 1:8)
cones <- FindClusters(cones, resolution = 2.0)
cones <- RunUMAP(object = cones, dims = 1:8,  reduction.name = "umap")


myRDS <- paste(mysample, "_SubsetCones.rds", sep="")
saveRDS(cones, file = myRDS)

figure_name <- ""
figure_name <- paste(mysample, "_ConesUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(cones, reduction = "umap",  repel = TRUE) + ggtitle("UMAP")
DimPlot(cones, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()

cell_values <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
mySubset <- subset(cones, idents = cell_values, invert = FALSE)

figure_name <- ""
figure_name <- paste(mysample, "_reClsusteredConesUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(mySubset, reduction = "umap",  repel = TRUE) + ggtitle("UMAP")
DimPlot(mySubset, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


mySubset@active.ident <- as.factor(mySubset@meta.data[,'sample'])
names(mySubset@active.ident) <- rownames(mySubset@meta.data)


myRDS <- paste(mysample, "_reClusteredCones.rds", sep="")
saveRDS(mySubset, file = myRDS)







